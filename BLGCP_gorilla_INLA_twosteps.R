rm(list = ls())
set.seed(2026)

library(INLA)
library(inlabru)
library(fmesher)
library(sf)
library(terra)
library(spatstat.geom)
library(spatstat.model)
library(spatstat.data)

data(gorillas_sf, package = "inlabru")
nests <- gorillas_sf$nests
mesh <- gorillas_sf$mesh
boundary <- gorillas_sf$boundary
nests_major <- nests[nests$group == "major", ]
nests_minor <- nests[nests$group == "minor", ]

data(gorillas, package = "spatstat.data")
X <- gorillas
W <- Window(X)
X_major <- unmark(X[X$marks$group == "major"])
X_minor <- unmark(X[X$marks$group == "minor"])

scale_im <- function(z) {
  stopifnot(inherits(z, "im"))
  vals <- as.vector(z$v)
  mu <- mean(vals, na.rm = TRUE)
  sdv <- sd(vals, na.rm = TRUE)
  if (!is.finite(sdv) || sdv <= 0) sdv <- 1
  z2 <- z
  z2$v <- (z$v - mu) / sdv
  attr(z2, "scaled:center") <- mu
  attr(z2, "scaled:scale") <- sdv
  z2
}

elevation_z_ppm  <- scale_im(gorillas.extra$elevation)
waterdist_z_ppm  <- scale_im(gorillas.extra$waterdist)
slopeangle_z_ppm <- scale_im(gorillas.extra$slopeangle)

cov_ppm <- list(
  elevation_z  = elevation_z_ppm,
  waterdist_z  = waterdist_z_ppm,
  slopeangle_z = slopeangle_z_ppm,
  heat         = gorillas.extra$heat,
  slopetype    = gorillas.extra$slopetype,
  vegetation   = gorillas.extra$vegetation
)

ppm_formula_major <- ~ elevation_z + waterdist_z + slopeangle_z +
  heat + vegetation

ppm_formula_minor <- ~ elevation_z + waterdist_z + slopeangle_z +
  heat + slopetype + vegetation

ppm_major <- ppm(
  X_major,
  trend = ppm_formula_major,
  covariates = cov_ppm
)

ppm_minor <- ppm(
  X_minor,
  trend = ppm_formula_minor,
  covariates = cov_ppm
)

beta_major <- coef(ppm_major)
beta_minor <- coef(ppm_minor)




im_to_spatraster_scaled <- function(z, scale_factor = 1000) {
  stopifnot(inherits(z, "im"))
  
  xcol <- z$xcol / scale_factor
  yrow <- z$yrow / scale_factor
  ny <- length(yrow)
  nx <- length(xcol)
  
  xyz <- data.frame(
    x = rep(xcol, each = ny),
    y = rep(yrow, times = nx),
    z = as.vector(z$v)   
  )
  
  terra::rast(xyz, type = "xyz")
}

trend_major_im <- predict(ppm_major, type = "trend")
trend_minor_im <- predict(ppm_minor, type = "trend")
logtrend_major_im <- trend_major_im
logtrend_minor_im <- trend_minor_im

logtrend_major_im$v <- log(pmax(trend_major_im$v, 1e-10))
logtrend_minor_im$v <- log(pmax(trend_minor_im$v, 1e-10))

logtrend_major <- im_to_spatraster_scaled(logtrend_major_im, scale_factor = 1000)
logtrend_minor <- im_to_spatraster_scaled(logtrend_minor_im, scale_factor = 1000)
gcov <- gorillas_sf_gcov()

terra::crs(logtrend_major) <- terra::crs(gcov$elevation)
terra::crs(logtrend_minor) <- terra::crs(gcov$elevation)
plot(is.na(logtrend_major))
plot(gorillas_sf$boundary, add = TRUE)




# exp SPDE models
# In 2D, alpha = nu + d/2. Exponential covariance means nu = 0.5, so alpha = 1.5.

spde_shared <- inla.spde2.pcmatern(
  mesh = gorillas_sf$mesh,
  alpha = 1.5,
  prior.range = c(0.1, 0.01),   
  prior.sigma = c(1, 0.01) 
)


spde_major <- inla.spde2.pcmatern(
  mesh = gorillas_sf$mesh,
  alpha = 1.5,
  prior.range = c(0.1, 0.01),
  prior.sigma = c(1, 0.01)
)

spde_minor <- inla.spde2.pcmatern(
  mesh = gorillas_sf$mesh,
  alpha = 1.5,
  prior.range = c(0.1, 0.01),
  prior.sigma = c(1, 0.01)
)


cmp <- ~
  Logtrend_major(logtrend_major, model = "offset") +
  Logtrend_minor(logtrend_minor, model = "offset") +
  Shared(geometry, model = spde_shared) +
  Major_field(geometry, model = spde_major) +
  Minor_field(geometry, model = spde_minor)

fml_major <- geometry ~ Logtrend_major + Shared + Major_field
fml_minor <- geometry ~ Logtrend_minor + Shared + Minor_field

lik_major <- bru_obs(
  "cp",
  formula = fml_major,
  data = nests_major,
  samplers = boundary,
  domain = list(geometry = mesh),
  tag = "major"
)

lik_minor <- bru_obs(
  "cp",
  formula = fml_minor,
  data = nests_minor,
  samplers = boundary,
  domain = list(geometry = mesh),
  tag = "minor"
)

t0 <- Sys.time()
fit_blgcp <- bru(
  cmp, lik_major, lik_minor,
  options = list(
    control.inla = list(
      #int.strategy = "eb"
      strategy = "auto",
      int.strategy = "auto"
    ),
    bru_max_iter = 10
  )
)
t1 <- Sys.time()
print(t1 - t0)


# Rerunning
fit_blgcp0 <- fit_blgcp
fit_blgcp <- bru_rerun(fit_blgcp)


print(fit_blgcp$summary.hyperpar)

# multiply by 1000 to return to the original coordinate scale. For alpha = 1.5/nu = 0.5 for exp,
# practical range is approximately 2 * xi.
# so x1000/2






