rm(list = ls())
set.seed(2026)

library(INLA)
library(inlabru)
library(fmesher)
library(sf)
library(terra)
data(gorillas_sf, package = "inlabru")
nests <- gorillas_sf$nests
mesh <- gorillas_sf$mesh
boundary <- gorillas_sf$boundary
gcov <- gorillas_sf_gcov()
nests_major <- nests[nests$group == "major", ]
nests_minor <- nests[nests$group == "minor", ]

standardise_spatraster <- function(r) {
  vals <- terra::values(r, mat = FALSE)
  mu <- mean(vals, na.rm = TRUE)
  sdv <- sd(vals, na.rm = TRUE)
  if (!is.finite(sdv) || sdv <= 0) sdv <- 1
  out <- (r - mu) / sdv
  list(raster = out, mean = mu, sd = sdv)
}

sc_elevation  <- standardise_spatraster(gcov$elevation)
sc_waterdist  <- standardise_spatraster(gcov$waterdist)
sc_slopeangle <- standardise_spatraster(gcov$slopeangle)
elevation_z  <- sc_elevation$raster
waterdist_z  <- sc_waterdist$raster
slopeangle_z <- sc_slopeangle$raster





spde_shared <- inla.spde2.pcmatern(
   mesh = gorillas_sf$mesh,
   alpha = 1.5,
   prior.range = c(0.1, 0.01),
   prior.sigma = c(1.0, 0.01)
 )
# 
spde_major <- inla.spde2.pcmatern(
  mesh = gorillas_sf$mesh,
  alpha = 1.5,
  prior.range = c(0.1, 0.01),
  prior.sigma = c(1.0, 0.01)
)
# 
spde_minor <- inla.spde2.pcmatern(
  mesh = gorillas_sf$mesh,
  alpha = 1.5,
  prior.range = c(0.1, 0.01),
  prior.sigma = c(1.0, 0.01)
)


cmp <- ~
  Shared(geometry, model = spde_shared) +
  Major_field(geometry, model = spde_major) +
  Minor_field(geometry, model = spde_minor) +
  Intercept_major(1) +
  Intercept_minor(1) +
  elevation_major(elevation_z, model = "linear") +
  waterdist_major(waterdist_z, model = "linear") +
  slopeangle_major(slopeangle_z, model = "linear") +
  elevation_minor(elevation_z, model = "linear") +
  waterdist_minor(waterdist_z, model = "linear") +
  slopeangle_minor(slopeangle_z, model = "linear") +
  heat_major(gcov$heat, model = "factor_contrast") +
  slopetype_major(gcov$slopetype, model = "factor_contrast") +
  vegetation_major(gcov$vegetation, model = "factor_contrast") +
  heat_minor(gcov$heat, model = "factor_contrast") +
  slopetype_minor(gcov$slopetype, model = "factor_contrast") +
  vegetation_minor(gcov$vegetation, model = "factor_contrast")


fml_major <- geometry ~ Intercept_major + elevation_major + waterdist_major + slopeangle_major + heat_major + slopetype_major + vegetation_major + Shared + Major_field
fml_minor <- geometry ~ Intercept_minor + elevation_minor + waterdist_minor + slopeangle_minor + heat_minor + slopetype_minor + vegetation_minor + Shared + Minor_field

lik_major <- bru_obs(
  "cp",
  formula = fml_major,
  data = gorillas_sf$nests[gorillas_sf$nests$group == "major", ],
  samplers = gorillas_sf$boundary,
  domain = list(geometry = gorillas_sf$mesh)
)

lik_minor <- bru_obs(
  "cp",
  formula = fml_minor,
  data = gorillas_sf$nests[gorillas_sf$nests$group == "minor", ],
  samplers = gorillas_sf$boundary,
  domain = list(geometry = gorillas_sf$mesh)
)

t0=Sys.time()
fit_blgcp <- bru(cmp, lik_major, lik_minor,
  options = list(
    control.inla = list(int.strategy = "eb"),
    bru_max_iter = 1
  )
)
t1=Sys.time()
t1-t0



# Rerunning
fit_blgcp0 <- fit_blgcp
fit_blgcp <- bru_rerun(fit_blgcp)



print(summary(fit_blgcp))
print(fit_blgcp$summary.fixed)
print(names(fit_blgcp$summary.random))
print(fit_blgcp$summary.hyperpar)



hyper <- fit_blgcp$summary.hyperpar
range_rows <- grep("^Range for", rownames(hyper))
xi_estimates <- hyper[range_rows, , drop = FALSE] / 2 *1000
rownames(xi_estimates) <- sub("Range for", "xi for", rownames(xi_estimates))
print(xi_estimates)





