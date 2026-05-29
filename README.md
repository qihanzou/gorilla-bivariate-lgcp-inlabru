# gorilla_bivariate_types_with_inlabru
Some inlabru code for gorilla data with two types of point patterns.

We assume a bivariate LGCP model:

$$
\log(\Lambda_i(s)) = \mu_i(s) + Y(s) + U_i(s),
$$

where $\mu_i(s)$ is the inhomogeneous mean trend with covariates, $Y(s)$ is the shared spatial field, and $U_i(s)$ is the individual spatial field for $i = 1, 2$.

1. `BLGCP_gorilla_INLA_full.R`: Full INLA/inlabru implementation for the bivariate LGCP model of the gorilla data with covariates.
2. `BLGCP_gorilla_INLA_twosteps.R`: Two-step method. The inhomogeneous mean trend is estimated by a Poisson model using `ppm()` in `spatstat`, and the latent fields are estimated using INLA/inlabru.
3. `INLA_env.R`: Functions for constructing the envelope test.
4. `INLA.R`: Additional utility functions.
