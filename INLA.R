
make_unit_square_boundary_sf = function(L = 1) {
  poly = matrix(
    c(0, 0, L, 0, L, L, 0, L, 0, 0),
    ncol = 2,
    byrow = TRUE
  )
  sf::st_sf(
    geometry = sf::st_sfc(sf::st_polygon(list(poly))),
    crs = NA
  )
}

ppp_to_sf_by_type = function(X, type_name) {
  Xi = X[marks(X) == type_name]
  
  if (Xi$n == 0) {
    return(sf::st_sf(
      geometry = sf::st_sfc(crs = NA)
    ))
  }
  
  df = data.frame(x = Xi$x, y = Xi$y)
  sf::st_as_sf(df, coords = c("x", "y"), crs = NA)
}

make_inla_prior_one = function(prior, name) {
  if (!is.null(prior[[name]])) {
    return(prior[[name]])
  }
  prior$default
}

fit_INLA_BLGCP_once = function(
    X,
    L = 1,
    mesh = NULL,
    boundary = NULL,
    prior = list(
      default = list(
        range = c(0.20, 0.50),
        sigma = c(1.00, 0.01)
      )
    ),
    mesh_max_edge = c(0.04, 0.12),
    mesh_cutoff = 0.005,
    int_strategy = "eb",
    bru_max_iter = 1,
    do_rerun = TRUE,
    verbose = FALSE
) {
  
  if (is.null(boundary)) {
    boundary = make_unit_square_boundary_sf(L)
  }
  
  if (is.null(mesh)) {
    mesh = fmesher::fm_mesh_2d(
      boundary = boundary,
      max.edge = mesh_max_edge,
      cutoff = mesh_cutoff
    )
  }
  
  dat_type1 = ppp_to_sf_by_type(X, "Type1")
  dat_type2 = ppp_to_sf_by_type(X, "Type2")
  
  if (nrow(dat_type1) == 0 || nrow(dat_type2) == 0) {
    return(list(ok = FALSE, est = rep(NA, 6), fit = NULL))
  }
  
  prior_shared = make_inla_prior_one(prior, "shared")
  prior_u1     = make_inla_prior_one(prior, "u1")
  prior_u2     = make_inla_prior_one(prior, "u2")
  
  spde_shared = INLA::inla.spde2.pcmatern(
    mesh = mesh,
    alpha = 1.5,
    prior.range = prior_shared$range,
    prior.sigma = prior_shared$sigma
  )
  
  spde_u1 = INLA::inla.spde2.pcmatern(
    mesh = mesh,
    alpha = 1.5,
    prior.range = prior_u1$range,
    prior.sigma = prior_u1$sigma
  )
  
  spde_u2 = INLA::inla.spde2.pcmatern(
    mesh = mesh,
    alpha = 1.5,
    prior.range = prior_u2$range,
    prior.sigma = prior_u2$sigma
  )
  
  cmp = ~
    Shared(geometry, model = spde_shared) +
    Type1_field(geometry, model = spde_u1) +
    Type2_field(geometry, model = spde_u2) +
    Intercept_type1(1) +
    Intercept_type2(1)
  
  fml_type1 = geometry ~ Intercept_type1 + Shared + Type1_field
  fml_type2 = geometry ~ Intercept_type2 + Shared + Type2_field
  
  lik_type1 = inlabru::bru_obs(
    "cp",
    formula = fml_type1,
    data = dat_type1,
    samplers = boundary,
    domain = list(geometry = mesh)
  )
  
  lik_type2 = inlabru::bru_obs(
    "cp",
    formula = fml_type2,
    data = dat_type2,
    samplers = boundary,
    domain = list(geometry = mesh)
  )
  
  fit = tryCatch(
    inlabru::bru(
      cmp,
      lik_type1,
      lik_type2,
      options = list(
        control.inla = list(int.strategy = int_strategy),
        bru_max_iter = bru_max_iter
      ),
      verbose = verbose
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(ok = FALSE, est = rep(NA, 6), fit = NULL))
  }
  
  if (isTRUE(do_rerun)) {
    fit = tryCatch(
      inlabru::bru_rerun(fit),
      error = function(e) fit
    )
  }
  
  hyper = fit$summary.hyperpar
  
  range_rows = grep("^Range for", rownames(hyper))
  sigma_rows = grep("^Stdev for", rownames(hyper))
  
  range_mean = hyper[range_rows, "mean"]
  sigma_mean = hyper[sigma_rows, "mean"]
  
  names(range_mean) = sub("^Range for ", "", names(range_mean))
  names(sigma_mean) = sub("^Stdev for ", "", names(sigma_mean))
  
  # INLA Matern nu = 0.5 range = 2 * exponential scale
  xi_shared = range_mean["Shared"] / 2
  xi_u1     = range_mean["Type1_field"] / 2
  xi_u2     = range_mean["Type2_field"] / 2
  
  sigma_shared = sigma_mean["Shared"]
  sigma_u1     = sigma_mean["Type1_field"]
  sigma_u2     = sigma_mean["Type2_field"]
  
  est = c(
    sigma_Y  = sigma_shared,
    scale_Y  = xi_shared,
    sigma_U1 = sigma_u1,
    sigma_U2 = sigma_u2,
    scale_U1 = xi_u1,
    scale_U2 = xi_u2
  )
  
  list(
    ok = all(is.finite(est)),
    est = est,
    fit = fit,
    hyper = hyper
  )
}