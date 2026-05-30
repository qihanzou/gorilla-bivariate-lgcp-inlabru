# Bivariate Gorilla Point Pattern Modeling with inlabru
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20457905.svg)](https://doi.org/10.5281/zenodo.20457905)

Some inlabru code for gorilla data with two types of point patterns.

We assume a bivariate LGCP model:

$$
\log(\Lambda_i(s)) = \mu_i(s) + Y(s) + U_i(s),
$$

where $\mu_i(s)$ is the inhomogeneous mean trend with covariates, $Y(s)$ is the shared spatial latent field, and $U_i(s)$ is the individual spatial latent field for $i = 1, 2$.

1. `BLGCP_gorilla_INLA_full.R`: Full INLA/inlabru implementation for the bivariate LGCP model of the gorilla data with covariates.
2. `BLGCP_gorilla_INLA_twosteps.R`: Two-step method. The inhomogeneous mean trend is estimated by a Poisson model using `ppm()` in `spatstat`, and the latent fields are estimated using INLA/inlabru.
3. `INLA_env.R`: Functions for constructing the envelope test.
4. `INLA.R`: Additional functions.
