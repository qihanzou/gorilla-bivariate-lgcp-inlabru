rm(list = ls())
set.seed(2026)
t0 = Sys.time()

library(spatstat.geom)
library(spatstat.explore)
library(spatstat.random)
library(spatstat.model)
library(spatstat.data)
library(reticulate)
library(future)
library(future.apply)
library(progressr)
library(abind)

image_nxy=128
data(gorillas, package = "spatstat.data")

X = gorillas
X_major = unmark(X[X$marks$group == "major"])
X_minor = unmark(X[X$marks$group == "minor"])
types = c("Type1", "Type2")

scale_im = function(z) {
  v = as.vector(z$v)
  eval.im((z - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
}

elevation = scale_im(gorillas.extra$elevation)
waterdist = scale_im(gorillas.extra$waterdist)
slopeangle = scale_im(gorillas.extra$slopeangle)
heat = gorillas.extra$heat
slopetype = gorillas.extra$slopetype
vegetation = gorillas.extra$vegetation

cov_list = list(
  elevation = elevation,
  waterdist = waterdist,
  slopeangle = slopeangle,
  heat = heat,
  slopetype = slopetype,
  vegetation = vegetation
)

fit_kppm_safe = function(Xi, cov_list, group = "major") {
  fit = tryCatch({
    if (group == "major") {
      kppm(unmark(Xi) ~ elevation + waterdist + slopeangle + heat + vegetation, clusters = "LGCP", data = cov_list)
    } else {
      kppm(unmark(Xi) ~ elevation + waterdist + slopeangle + heat + slopetype + vegetation, clusters = "LGCP", data = cov_list)
    }
  }, error = function(e) NULL)
  fit
}




get_lambda_trend = function(fit, beta = NULL, dimyx = c(50, 50)) {
  fit_ppm = as.ppm(fit)
  if (is.null(beta)) {
    out = predict(fit_ppm, type = "trend", dimyx = dimyx)
  } else {
    names(beta) = names(coef(fit_ppm))
    out = predict(fit_ppm, type = "trend", new.coef = beta, dimyx = dimyx)
  }
  out
}

get_beta_vec = function(fit1, fit2) {
  b1 = coef(fit1)
  b2 = coef(fit2)
  c(
    b1["(Intercept)"], b2["(Intercept)"],
    b1["elevation"], b2["elevation"],
    b1["waterdist"], b2["waterdist"],
    b1["slopeangle"], b2["slopeangle"],
    b1["heatModerate"], b2["heatModerate"],
    b1["heatCoolest"], b2["heatCoolest"],
    b2["slopetypeToe"],
    b2["slopetypeFlat"],
    b2["slopetypeMidslope"],
    b2["slopetypeUpper"],
    b2["slopetypeRidge"],
    b1["vegetationColonising"], b2["vegetationColonising"],
    b1["vegetationGrassland"], b2["vegetationGrassland"],
    b1["vegetationPrimary"], b2["vegetationPrimary"],
    b1["vegetationSecondary"], b2["vegetationSecondary"],
    b1["vegetationTransition"], b2["vegetationTransition"]
  )
}


scale_by_train = function(train_x, test_x) {
  scale_mean = colMeans(train_x, na.rm = TRUE)
  scale_sd = apply(train_x, 2, sd, na.rm = TRUE)
  scale_sd[scale_sd == 0] = 1
  train_scaled = sweep(train_x, 2, scale_mean, "-")
  train_scaled = sweep(train_scaled, 2, scale_sd, "/")
  test_scaled = sweep(test_x, 2, scale_mean, "-")
  test_scaled = sweep(test_scaled, 2, scale_sd, "/")
  list(train = train_scaled, test = test_scaled, mean = scale_mean, sd = scale_sd)
}

fit_major = fit_kppm_safe(X_major, cov_list, "major")
fit_minor = fit_kppm_safe(X_minor, cov_list, "minor")
ppm_major = as.ppm(fit_major)
ppm_minor = as.ppm(fit_minor)
beta_major = coef(ppm_major)
beta_minor = coef(ppm_minor)


lambda_mean_trend_major = get_lambda_trend(fit_major, dimyx = c(image_nxy, image_nxy))
lambda_mean_trend_minor = get_lambda_trend(fit_minor, dimyx = c(image_nxy, image_nxy))
mean_trend_major = eval.im(log(lambda_mean_trend_major))
mean_trend_minor = eval.im(log(lambda_mean_trend_minor))
W_major = Window(mean_trend_major)
W_minor = Window(mean_trend_minor)
W = union.owin(W_major, W_minor)
zero_mu = eval.im(0 * mean_trend_major)
X_obs = superimpose(Type1 = X_major, Type2 = X_minor, W = W)
marks(X_obs) = factor(marks(X_obs), levels = types)
rmax = rmax.rule("K", W)
r = seq(0, rmax, length.out = 513)





env_mlgcp_cross_fast <- function(theta, method_name = "DSBI",
                                 nsim = 99, r_use = r,
                                 max_tries = 100) {
  sigma_y  <- theta[1]
  scale_y  <- theta[2]
  sigma_u1 <- theta[3]
  sigma_u2 <- theta[4]
  scale_u1 <- theta[5]
  scale_u2 <- theta[6]
  
  env_cross <- envelope(
    X_obs,
    fun = function(Y, r) {
      Kcross.inhom(
        Y, i = "Type1", j = "Type2", r = r,
        lambdaI = lambda_mean_trend_major,
        lambdaJ = lambda_mean_trend_minor,
        correction = "border"
      )
    },
    r = r_use,
    simulate = expression({
      Y_shared <- log(attr(rLGCP(
        model = "exponential",
        mu = zero_mu,
        var = sigma_y^2,
        scale = scale_y,
        win = W,
        saveLambda = TRUE
      ), "Lambda"))
      
      mu_major <- eval.im(
        mean_trend_major -
          0.5 * sigma_y^2 -
          0.5 * sigma_u1^2 +
          Y_shared
      )
      
      mu_minor <- eval.im(
        mean_trend_minor -
          0.5 * sigma_y^2 -
          0.5 * sigma_u2^2 +
          Y_shared
      )
      
      ok <- FALSE
      tries <- 0
      
      while (!ok && tries < max_tries) {
        tries <- tries + 1
        
        X_major_sim <- rLGCP(
          model = "exponential",
          mu = mu_major,
          win = W,
          var = sigma_u1^2,
          scale = scale_u1,
          saveLambda = FALSE
        )
        
        X_minor_sim <- rLGCP(
          model = "exponential",
          mu = mu_minor,
          win = W,
          var = sigma_u2^2,
          scale = scale_u2,
          saveLambda = FALSE
        )
        
        ok <- npoints(X_major_sim) > 0 && npoints(X_minor_sim) > 0
      }
      
      if (!ok) {
        stop("Failed to simulate both Type1 and Type2 points after max_tries.")
      }
      
      Xsim <- superimpose(Type1 = X_major_sim, Type2 = X_minor_sim, W = W)
      marks(Xsim) <- factor(marks(Xsim), levels = c("Type1", "Type2"))
      Xsim
    }),
    nsim = nsim,
    savefuns = TRUE,
    global = FALSE,
    verbose = FALSE
  )
  
  test_cross <- dclf.test(env_cross)
  plot(env_cross, main = paste0("Cross-type envelope (", method_name, ")"))
  
  list(env = env_cross, test = test_cross)
}


env_mlgcp_major_fast <- function(theta, method_name = "DSBI",
                                 nsim = 99, r_use = r,
                                 max_tries = 100) {
  sigma_y  <- theta[1]
  scale_y  <- theta[2]
  sigma_u1 <- theta[3]
  scale_u1 <- theta[5]
  
  X_obs_major <- unmark(X_obs[marks(X_obs) == "Type1"])
  
  env_major <- envelope(
    X_obs_major,
    fun = function(Y, r) {
      Kinhom(
        Y, r = r,
        lambda = lambda_mean_trend_major,
        correction = "border"
      )
    },
    r = r_use,
    simulate = expression({
      Y_shared <- log(attr(rLGCP(
        model = "exponential",
        mu = zero_mu,
        var = sigma_y^2,
        scale = scale_y,
        win = W,
        saveLambda = TRUE
      ), "Lambda"))
      
      mu_major <- eval.im(
        mean_trend_major -
          0.5 * sigma_y^2 -
          0.5 * sigma_u1^2 +
          Y_shared
      )
      
      ok <- FALSE
      tries <- 0
      
      while (!ok && tries < max_tries) {
        tries <- tries + 1
        
        X_major_sim <- rLGCP(
          model = "exponential",
          mu = mu_major,
          win = W,
          var = sigma_u1^2,
          scale = scale_u1,
          saveLambda = FALSE
        )
        
        ok <- npoints(X_major_sim) > 0
      }
      
      if (!ok) {
        stop("Failed to simulate Major points after max_tries.")
      }
      
      X_major_sim
    }),
    nsim = nsim,
    savefuns = TRUE,
    global = FALSE,
    verbose = FALSE
  )
  
  test_major <- dclf.test(env_major)
  plot(env_major, main = paste0("Major group envelope (", method_name, ")"))
  
  list(env = env_major, test = test_major)
}



env_mlgcp_minor_fast <- function(theta, method_name = "DSBI",
                                 nsim = 99, r_use = r,
                                 max_tries = 100) {
  sigma_y  <- theta[1]
  scale_y  <- theta[2]
  sigma_u2 <- theta[4]
  scale_u2 <- theta[6]
  
  X_obs_minor <- unmark(X_obs[marks(X_obs) == "Type2"])
  
  env_minor <- envelope(
    X_obs_minor,
    fun = function(Y, r) {
      Kinhom(
        Y, r = r,
        lambda = lambda_mean_trend_minor,
        correction = "border"
      )
    },
    r = r_use,
    simulate = expression({
      Y_shared <- log(attr(rLGCP(
        model = "exponential",
        mu = zero_mu,
        var = sigma_y^2,
        scale = scale_y,
        win = W,
        saveLambda = TRUE
      ), "Lambda"))
      
      mu_minor <- eval.im(
        mean_trend_minor -
          0.5 * sigma_y^2 -
          0.5 * sigma_u2^2 +
          Y_shared
      )
      
      ok <- FALSE
      tries <- 0
      
      while (!ok && tries < max_tries) {
        tries <- tries + 1
        
        X_minor_sim <- rLGCP(
          model = "exponential",
          mu = mu_minor,
          win = W,
          var = sigma_u2^2,
          scale = scale_u2,
          saveLambda = FALSE
        )
        
        ok <- npoints(X_minor_sim) > 0
      }
      
      if (!ok) {
        stop("Failed to simulate Minor points after max_tries.")
      }
      
      X_minor_sim
    }),
    nsim = nsim,
    savefuns = TRUE,
    global = FALSE,
    verbose = FALSE
  )
  
  test_minor <- dclf.test(env_minor)
  plot(env_minor, main = paste0("Minor group envelope (", method_name, ")"))
  
  list(env = env_minor, test = test_minor)
}












theta_hat = c(sigma_y, scale_y, sigma_u1, sigma_u2, scale_u1, scale_u2)



set.seed(2025)
plot.res = 600
png("gorilla_env_cross_DSBI.png", width = 6.2 * plot.res, height = 4.5 * plot.res, res = plot.res)
res_cross = env_mlgcp_cross_fast(theta_hat, "DSBI")
dev.off()
set.seed(2025)
png("gorilla_env_major_DSBI.png", width = 6.2 * plot.res, height = 4.5 * plot.res, res = plot.res)
res_major = env_mlgcp_major_fast(theta_hat, "DSBI")
dev.off()
set.seed(2025)
png("gorilla_env_minor_DSBI.png", width = 6.2 * plot.res, height = 4.5 * plot.res, res = plot.res)
res_minor = env_mlgcp_minor_fast(theta_hat, "DSBI")
dev.off()

save.image("gorilla_DSBI_application_workspace.RData")
Sys.time() - t0



res_cross$env
res_major$env
res_minor$env

library(GET)
res_get1 <- global_envelope_test(curve_sets = list(Cross = res_cross$env, Major = res_major$env, Minor = res_minor$env), type = "erl")
res_get1
plot(res_get1)


plot.res <- 600
png("G_env.png", width = 6.2*plot.res, height = 4.5*plot.res, res = plot.res)
p <- plot(res_get1)
dev.off()





