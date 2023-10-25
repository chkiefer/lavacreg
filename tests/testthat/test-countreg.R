# ---------------------------------------------------
# TEST 1 - intercept only - one group
# ---------------------------------------------------
test_that("intercept-only Poisson", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ 1",
    group = NULL,
    data = example01,
    family = "poisson"
  )
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 3)
  expect_equal(length(avar), 2)

  # LOG-LIKELIHOOD
  comp <- 3.38334
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- 6.76966
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5)

  # 2. Regression coefficient
  comp <- 2.55546
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5)

  # 3. Overdispersion parameter
  comp <- 0
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5)

  # STANDARD ERRORS
  # 1. Group weight
  comp <- 0.00115
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5)

  # 2. Regression coefficient
  comp <- 0.00009
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5)
})


test_that("intercept-only negative binomial", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ 1",
    group = NULL,
    data = example01,
    family = "nbinom"
  )
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 3)
  expect_equal(length(avar), 3)

  # LOG-LIKELIHOOD
  comp <- 3.12601
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- 6.76961
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  comp <- 2.55545
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 3. Overdispersion parameter
  comp <- 8.62666
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # STANDARD ERRORS
  # 1. Group weight
  comp <- 0.00115
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  comp <- 0.00022
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 3. Overdispersion parameter
  comp <- 0.51543
  par <- avar[pt$par_free[pt$dest == "overdis"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")
})


# ---------------------------------------------------
# TEST 2 - intercept-only - two groups
# ---------------------------------------------------
test_that("two-group intercept-only Poisson", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ 1",
    group = "treat",
    data = example01,
    family = "poisson"
  )
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 6)
  expect_equal(length(avar), 4)

  # LOG-LIKELIHOOD
  comp <- 3.31971
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  comp <- c(2.44539, 2.65140)
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 3. Overdispersion parameter
  comp <- c(0, 0)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  comp <- c(0.00020, 0.00016)
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")
})


test_that("two-group intercept-only negative binomial", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ 1",
    group = "treat",
    data = example01,
    family = "nbinom"
  )
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 6)
  expect_equal(length(avar), 6)

  # LOG-LIKELIHOOD
  comp <- 3.09931
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  comp <- c(2.44539, 2.65140)
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 3. Overdispersion parameter
  comp <- c(7.90338, 11.58069)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00233, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  comp <- c(0.00050, 0.00035)
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 3. Overdispersion parameter
  comp <- c(0.88228, 2.20150)
  par <- avar[pt$par_free[pt$dest == "overdis"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")
})


# ---------------------------------------------------
# TEST 3 - one manifest covariate - two groups
# ---------------------------------------------------
test_that("two-group one manifest covariate Poisson", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ z12",
    group = "treat",
    data = example01,
    family = "poisson",
    se = TRUE
  )
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 12)
  expect_equal(length(avar), 10)

  # LOG-LIKELIHOOD
  comp <- 4.73872
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  comp <- c(2.65523, -0.16981, 2.83597, -0.14146)
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 3. Overdispersion parameter
  comp <- c(0, 0)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(1.36293, 1.40557)
  par <- pt$par[pt$dest == "mu_z"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

  # 4.2 Variances
  comp <- c(1.59051, 1.48010)
  par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")


  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  comp <- c(0.00039, 0.00015, 0.00034, 0.00012)
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(0.00372, 0.00334)
  par <- avar[pt$par_free[pt$dest == "mu_z"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

  # 4.2 Variances
  comp <- c(0.01182, 0.00989)
  par <- avar[pt$par_free[pt$dest == "Sigma_z"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")
})


test_that("two-group one manifest covariate negative binomial", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ z12",
    group = "treat",
    data = example01,
    family = "nbinom",
    se = TRUE
  )
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 12)
  expect_equal(length(avar), 12)

  # LOG-LIKELIHOOD
  comp <- 4.62806
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  comp <- c(2.65465, -0.16934, 2.83624, -0.14166)
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 3. Overdispersion parameter
  comp <- c(12.08287, 17.89945)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(1.36293, 1.40557)
  par <- pt$par[pt$dest == "mu_z"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

  # 4.2 Variances
  comp <- c(1.59050, 1.48011)
  par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")


  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  comp <- c(0.00082, 0.00028, 0.00063, 0.00021)
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 3. Overdispersion parameter
  comp <- c(3.12545, 8.45198)
  par <- avar[pt$par_free[pt$dest == "overdis"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(0.00372, 0.00334)
  par <- avar[pt$par_free[pt$dest == "mu_z"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

  # 4.2 Variances
  comp <- c(0.01182, 0.00989)
  par <- avar[pt$par_free[pt$dest == "Sigma_z"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")
})


# ---------------------------------------------------
# TEST 4 - three manifest covariates - two groups
# ---------------------------------------------------
test_that("two-group three manifest covariates Poisson", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ z12 + z11 + z21",
    group = "treat",
    data = example01,
    family = "poisson",
    se = TRUE
  )
  # Converged?
  conv <- fit@fit@fit$convergence
  expect_equal(conv, 0)

  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 30)
  expect_equal(length(avar), 28)

  # LOG-LIKELIHOOD
  comp <- 7.27241
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  comp <- c(
    2.34110, -0.08721, -0.07893, 0.08227,
    2.62149, -0.09299, -0.04680, 0.05441
  )
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 3. Overdispersion parameter
  comp <- c(0, 0)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(
    1.36293, 1.59813, 3.90927,
    1.40558, 1.56209, 3.99473
  )
  par <- pt$par[pt$dest == "mu_z"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

  # 4.2 Variances
  comp <- c(
    1.59049, 1.68417, 0.95031,
    1.48009, 1.36278, 0.83293
  )
  par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

  # 4.3 Covariances
  comp <- c(
    1.02240, -0.63361, -0.50519,
    0.86383, -0.56871, -0.46656
  )
  par <- pt$par[pt$dest == "Sigma_z" & pt$type == "cov"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_cov")


  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  comp <- c(
    0.00763, 0.00027, 0.00021, 0.00034,
    0.00711, 0.00021, 0.00020, 0.00030
  )
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(
    0.00372, 0.00394, 0.00222,
    0.00334, 0.00307, 0.00188
  )
  par <- avar[pt$par_free[pt$dest == "mu_z"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

  # 4.2 Variances
  comp <- c(
    0.01182, 0.01325, 0.00422,
    0.00989, 0.00838, 0.00313
  )
  par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")

  # 4.3 Covariances
  comp <- c(
    0.00870, 0.00447, 0.00434,
    0.00624, 0.00351, 0.00305
  )
  par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "cov"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_cov")
})


test_that("two-group three manifest covariates negative binomial", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ z12 + z11 + z21",
    group = "treat",
    data = example01,
    family = "nbinom",
    se = TRUE
  )
  # Converged?
  conv <- fit@fit@fit$convergence
  expect_equal(conv, 0)

  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 30)
  expect_equal(length(avar), 30)

  # LOG-LIKELIHOOD
  comp <- 7.18300
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  comp <- c(
    2.33834, -0.08674, -0.07975, 0.08312,
    2.63161, -0.09387, -0.04665, 0.05214
  )
  par <- pt$par[pt$dest == "beta"] |> round(5)
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 3. Overdispersion parameter
  comp <- c(14.17296, 19.56794)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(
    1.36293, 1.59813, 3.90927,
    1.40557, 1.56208, 3.99473
  )
  par <- pt$par[pt$dest == "mu_z"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

  # 4.2 Variances
  comp <- c(
    1.59051, 1.68417, 0.95030,
    1.48011, 1.36278, 0.83293
  )
  par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

  # 4.3 Covariances
  comp <- c(
    1.02241, -0.63361, -0.50518,
    0.86384, -0.56872, -0.46656
  )
  par <- pt$par[pt$dest == "Sigma_z" & pt$type == "cov"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_cov")


  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  comp <- c(
    0.01303, 0.00047, 0.00037, 0.00058,
    0.01184, 0.00036, 0.00034, 0.00050
  )
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 3. Overdispersion parameter
  comp <- c(5.07994, 11.22616)
  par <- avar[pt$par_free[pt$dest == "overdis"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(
    0.00372, 0.00393, 0.00222,
    0.00334, 0.00308, 0.00188
  )
  par <- avar[pt$par_free[pt$dest == "mu_z"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

  # 4.2 Variances
  comp <- c(
    0.01182, 0.01325, 0.00422,
    0.00989, 0.00838, 0.00313
  )
  par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")

  # 4.3 Covariances
  comp <- c(
    0.00870, 0.00447, 0.00434,
    0.00624, 0.00351, 0.00305
  )
  par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "cov"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_cov")
})


# ---------------------------------------------------
# TEST 5 - one latent covariate - two groups
# ---------------------------------------------------
test_that("two-group one latent covariate Poisson", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ eta1",
    lv = list(eta1 = c("z21", "z22")),
    group = "treat",
    data = example01,
    family = "poisson",
    se = TRUE
  )
  # Converged?
  conv <- fit@fit@fit$convergence
  expect_equal(conv, 0L)

  # Correct parameter estimates?
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 24)
  expect_equal(length(avar), 18)

  # LOG-LIKELIHOOD
  comp <- 5.41229
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(0.42565, 0.81313)
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 2.2 gamma
  comp <- c(0.50043, 0.45471)
  par <- pt$par[pt$dest == "gamma"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_gamma")

  # 3. Overdispersion parameter
  comp <- c(0, 0)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0, 1.45445), 2)
  par <- pt$par[pt$dest == "nu"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(1, 0.72492), 2)
  par <- pt$par[pt$dest == "Lambda"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_lambda")

  # 5.3 Theta
  comp <- c(0.52596, 0.35650, 0.50998, 0.35629)
  par <- pt$par[pt$dest == "Theta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(3.93444, 3.97110)
  par <- pt$par[pt$dest == "mu_eta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_eta")

  # 6.2 Variances
  comp <- c(0.42556, 0.32470)
  par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_var")

  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(0.04255, 0.04600)
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 2.2 gamma
  comp <- c(0.00258, 0.00281)
  par <- avar[pt$par_free[pt$dest == "gamma"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_gamma")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(0.04933, 2)
  par <- avar[pt$par_free[pt$dest == "nu"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_nu")

  # 5.2 Lambda
  comp <- rep(0.00311, 2)
  par <- avar[pt$par_free[pt$dest == "Lambda"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_lambda")

  # 5.3 Theta
  comp <- c(0.00292, 0.00111, 0.00261, 0.00101)
  par <- avar[pt$par_free[pt$dest == "Theta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(0.00196, 0.00165)
  par <- avar[pt$par_free[pt$dest == "mu_eta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_eta")

  # 6.2 Variances
  comp <- c(0.00404, 0.00285)
  par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_var")
})


test_that("two-group one latent covariate negative binomial", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ eta1",
    lv = list(eta1 = c("z21", "z22")),
    group = "treat",
    data = example01,
    family = "nbinom",
    se = TRUE
  )
  # Converged?
  conv <- fit@fit@fit$convergence
  expect_equal(conv, 0)

  # Correct parameter estimates?
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 24)
  expect_equal(length(avar), 20)

  # LOG-LIKELIHOOD
  comp <- 5.36203
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(1.39488, 1.38017)
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 2.2 gamma
  comp <- c(0.26256, 0.31877)
  par <- pt$par[pt$dest == "gamma"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_gamma")

  # 3. Overdispersion parameter
  comp <- c(9.58424, 15.98937)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0, -0.74553), 2)
  par <- pt$par[pt$dest == "nu"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(1, 1.2814), 2)
  par <- pt$par[pt$dest == "Lambda"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_lambda")

  # 5.3 Theta
  comp <- c(0.66769, 0.03286, 0.57080, 0.16991)
  par <- pt$par[pt$dest == "Theta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(3.87379, 3.95331)
  par <- pt$par[pt$dest == "mu_eta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_eta")

  # 6.2 Variances
  comp <- c(0.45845, 0.22501)
  par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_var")


  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(0.03008, 0.04599)
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 2.2 gamma
  comp <- c(0.00185, 0.00285)
  par <- avar[pt$par_free[pt$dest == "gamma"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_gamma")

  # 3. Overdispersion parameter
  comp <- c(1.57327, 7.36865)
  par <- avar[pt$par_free[pt$dest == "overdis"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(0.12104, 2)
  par <- avar[pt$par_free[pt$dest == "nu"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_nu")

  # 5.2 Lambda
  comp <- rep(0.00765, 2)
  par <- avar[pt$par_free[pt$dest == "Lambda"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_lambda")

  # 5.3 Theta
  comp <- c(0.00215, 0.00001, 0.00221, 0.00181)
  par <- avar[pt$par_free[pt$dest == "Theta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(0.00088, 0.00127)
  par <- avar[pt$par_free[pt$dest == "mu_eta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_eta")

  # 6.2 Variances
  comp <- c(0.00389, 0.00098)
  par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_var")
})


# ---------------------------------------------------
# TEST 6 - two latent covariates - two groups
# ---------------------------------------------------
test_that("two-group two latent covariates Poisson", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ eta1 + eta2",
    lv = list(
      eta1 = c("z21", "z22"),
      eta2 = c("z41", "z42", "z43")
    ),
    group = "treat",
    data = example01,
    family = "poisson",
    se = TRUE
  )
  # Converged?
  conv <- fit@fit@fit$convergence
  expect_equal(conv, 0)

  # Correct parameter estimates?
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 60)
  expect_equal(length(avar), 40)

  # LOG-LIKELIHOOD
  comp <- 11.09449
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")


  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(0.63028, 0.55231)
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 2.2 gamma
  comp <- c(0.46326, -0.03616, 0.50796, 0.02916)
  par <- pt$par[pt$dest == "gamma"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_gamma")

  # 3. Overdispersion parameter
  comp <- c(0, 0)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0, 1.35728, 0, -0.09500, -0.45658), 2)
  par <- pt$par[pt$dest == "nu"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(1, 0.74971, rep(0, 5), 1, 1.28107, 1.36110), 2)
  par <- pt$par[pt$dest == "Lambda"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_lambda")

  # 5.3 Theta
  comp <- c(
    0.52566, 0.33860, 1.53647, 1.43656, 1.39938,
    0.52605, 0.35564, 1.36698, 1.56639, 1.02974
  )
  par <- pt$par[pt$dest == "Theta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(3.93991, 1.57989, 3.97623, 1.65178)
  par <- pt$par[pt$dest == "mu_eta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_eta")

  # 6.2 Variances
  comp <- c(0.43188, 1.89513, 0.31349, 2.16220)
  par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_var")

  # 6.3 Covariances
  comp <- c(-0.34766, -0.41424)
  par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "cov"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_cov")


  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(0.06556, 0.09605)
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 2.2 gamma
  comp <- c(0.00340, 0.00044, 0.00508, 0.00040)
  par <- avar[pt$par_free[pt$dest == "gamma"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_gamma")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0.04997, 0.01275, 0.01540), 2)
  par <- avar[pt$par_free[pt$dest == "nu"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(0.00315, 0.00307, 0.00404), 2)
  par <- avar[pt$par_free[pt$dest == "Lambda"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_lambda")

  # 5.3 Theta
  comp <- c(
    0.00314, 0.00109, 0.01784, 0.02554, 0.02826,
    0.00256, 0.00099, 0.01433, 0.02544, 0.03841
  )
  par <- avar[pt$par_free[pt$dest == "Theta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(0.00194, 0.00606, 0.00156, 0.00518)
  par <- avar[pt$par_free[pt$dest == "mu_eta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_eta")

  # 6.2 Variances
  comp <- c(0.00431, 0.04134, 0.00277, 0.05227)
  par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_var")

  # 6.3 Covariances
  comp <- c(0.00478, 0.00470)
  par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "cov"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_cov")
})



test_that("two-group two latent covariates negative binomial", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ eta1 + eta2",
    lv = list(
      eta1 = c("z21", "z22"),
      eta2 = c("z41", "z42", "z43")
    ),
    group = "treat",
    data = example01,
    family = "nbinom",
    se = TRUE
  )
  # Converged?
  conv <- fit@fit@fit$convergence
  expect_equal(conv, 0)

  # Correct parameter estimates?
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 60)
  expect_equal(length(avar), 42)

  # LOG-LIKELIHOOD
  comp <- 11.03596
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(1.42405, 1.53285)
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 2.2 gamma
  comp <- c(0.28020, -0.06825, 0.28419, -0.01419)
  par <- pt$par[pt$dest == "gamma"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_gamma")

  # 3. Overdispersion parameter
  comp <- c(14.00803, 18.66380)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0, 1.08882, 0, -0.07963, -0.73781), 2)
  par <- pt$par[pt$dest == "nu"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(1, 0.81774, rep(0, 5), 1, 1.27078, 1.53127), 2)
  par <- pt$par[pt$dest == "Lambda"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_lambda")

  # 5.3 Theta
  comp <- c(
    0.45836, 0.25828, 1.61553, 1.61485, 1.17764,
    0.45605, 0.28609, 1.68552, 2.03674, 0.28088
  )
  par <- pt$par[pt$dest == "Theta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(3.94419, 1.58128, 3.96585, 1.684321)
  par <- pt$par[pt$dest == "mu_eta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_eta")

  # 6.2 Variances
  comp <- c(0.48812, 1.66195, 0.34983, 2.08168)
  par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_var")

  # 6.3 Covariances
  comp <- c(-0.30853, -0.30940)
  par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "cov"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_cov")


  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(0.03527, 0.04232)
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 2.2 gamma
  comp <- c(0.00184, 0.00041, 0.00227, 0.00024)
  par <- avar[pt$par_free[pt$dest == "gamma"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_gamma")

  # 3. Overdispersion parameter
  comp <- c(6.3742, 12.7716)
  par <- avar[pt$par_free[pt$dest == "overdis"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0.09798, 0.01402, 0.01389), 2)
  par <- avar[pt$par_free[pt$dest == "nu"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(0.00622, 0.00330, 0.00322), 2)
  par <- avar[pt$par_free[pt$dest == "Lambda"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_lambda")

  # 5.3 Theta
  comp <- c(
    0.00353, 0.00140, 0.01872, 0.02576, 0.03203,
    0.00250, 0.00110, 0.01433, 0.02167, 0.00095
  )
  par <- avar[pt$par_free[pt$dest == "Theta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(0.00192, 0.00559, 0.00145, 0.00276)
  par <- avar[pt$par_free[pt$dest == "mu_eta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_eta")

  # 6.2 Variances
  comp <- c(0.00496, 0.02802, 0.00237, 0.02454)
  par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_var")

  # 6.3 Covariances
  comp <- c(0.00399, 0.00287)
  par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "cov"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_cov")
})


# ---------------------------------------------------
# TEST 7 - one latent, one manifest covariates - two groups
# ---------------------------------------------------
test_that("two-group one latent, one manifest covariate Poisson", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ eta1 + z12",
    lv = list(eta1 = c("z41", "z42", "z43")),
    group = "treat",
    data = example01,
    family = "poisson",
    se = TRUE
  )

  # Converged?
  conv <- fit@fit@fit$convergence
  expect_equal(conv, 0)

  # Correct parameter estimates?
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 38)
  expect_equal(length(avar), 32)

  # LOG-LIKELIHOOD
  comp <- 10.40298
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")


  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(2.76045, -0.14176, 2.86583, -0.12890)
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 2.2 gamma
  comp <- c(-0.09412, -0.02862)
  par <- pt$par[pt$dest == "gamma"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_gamma")

  # 3. Overdispersion parameter
  comp <- c(0, 0)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(1.35905, 1.39256)
  par <- pt$par[pt$dest == "mu_z"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

  # 4.2 Variances
  comp <- c(1.59276, 1.49301)
  par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0, -0.09046, -0.43259), 2)
  par <- pt$par[pt$dest == "nu"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(1, 1.27837, 1.34651), 2)
  par <- pt$par[pt$dest == "Lambda"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_lambda")

  # 5.3 Theta
  comp <- c(
    1.51950, 1.46891, 1.45808,
    1.36367, 1.51842, 1.08572
  )
  par <- pt$par[pt$dest == "Theta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(1.58518, 1.64336)
  par <- pt$par[pt$dest == "mu_eta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_eta")

  # 6.2 Variances
  comp <- c(1.87179, 2.13073)
  par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_var")

  # 7. Latent-Manifest Covariances
  comp <- c(0.50240, 0.68700)
  par <- pt$par[pt$dest == "Sigma_z_lv"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_eta_cov")


  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(0.00070, 0.00018, 0.00048, 0.00015)
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 2.2 gamma
  comp <- c(0.00024, 0.00014)
  par <- avar[pt$par_free[pt$dest == "gamma"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_gamma")

  # # 3. Overdispersion parameter
  # comp <- c(9.21046, 7.77312)
  # par <- avar[pt$par_free[pt$dest == "overdis"]]
  # expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(0.00370, 0.00323)
  par <- avar[pt$par_free[pt$dest == "mu_z"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

  # 4.2 Variances
  comp <- c(0.01190, 0.01034)
  par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")


  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0.01280, 0.01482), 2)
  par <- avar[pt$par_free[pt$dest == "nu"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(0.00310, 0.00384), 2)
  par <- avar[pt$par_free[pt$dest == "Lambda"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_lambda")

  # 5.3 Theta
  comp <- c(
    0.01775, 0.02544, 0.02820,
    0.01410, 0.02401, 0.03513
  )
  par <- avar[pt$par_free[pt$dest == "Theta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(0.00625, 0.00539)
  par <- avar[pt$par_free[pt$dest == "mu_eta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_eta")

  # 6.2 Variances
  comp <- c(0.03885, 0.05531)
  par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_var")

  # 7. Latent-manifest Covariances
  comp <- c(0.00957, 0.01124)
  par <- avar[pt$par_free[pt$dest == "Sigma_z_lv"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_cov")
})


test_that("two-group one latent, one manifest covariate negative binomial", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ eta1 + z12",
    lv = list(eta1 = c("z41", "z42", "z43")),
    group = "treat",
    data = example01,
    family = "nbinom",
    se = TRUE
  )
  # Converged?
  conv <- fit@fit@fit$convergence
  expect_equal(conv, 0)

  # Correct parameter estimates?
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 38)
  expect_equal(length(avar), 34)

  # LOG-LIKELIHOOD
  comp <- 10.30759
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(2.74023, -0.14692, 2.86376, -0.13065)
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 2.2 gamma
  comp <- c(-0.07545, -0.02587)
  par <- pt$par[pt$dest == "gamma"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_gamma")

  # 3. Overdispersion parameter
  comp <- c(13.57181, 18.25335)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(1.35814, 1.39329)
  par <- pt$par[pt$dest == "mu_z"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

  # 4.2 Variances
  comp <- c(1.59198, 1.49435)
  par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0, -0.10521, -0.46664), 2)
  par <- pt$par[pt$dest == "nu"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(1, 1.28735, 1.36718), 2)
  par <- pt$par[pt$dest == "Lambda"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_lambda")

  # 5.3 Theta
  comp <- c(
    1.55109, 1.45711, 1.36321,
    1.38212, 1.52278, 1.04421
  )
  par <- pt$par[pt$dest == "Theta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(1.58199, 1.64565)
  par <- pt$par[pt$dest == "mu_eta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_eta")

  # 6.2 Variances
  comp <- c(1.84814, 2.11163)
  par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_var")

  # 7. Latent-Manifest Covariances
  comp <- c(0.49467, 0.68556)
  par <- pt$par[pt$dest == "Sigma_z_lv"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_eta_cov")


  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(0.00119, 0.00030, 0.00086, 0.00024)
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 2.2 gamma
  comp <- c(0.00031, 0.00020)
  par <- avar[pt$par_free[pt$dest == "gamma"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_gamma")

  # 3. Overdispersion parameter
  comp <- c(4.58670, 8.98671)
  par <- avar[pt$par_free[pt$dest == "overdis"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(0.00370, 0.00322)
  par <- avar[pt$par_free[pt$dest == "mu_z"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

  # 4.2 Variances
  comp <- c(0.01187, 0.01038)
  par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")


  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0.01298, 0.01573), 2)
  par <- avar[pt$par_free[pt$dest == "nu"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(0.00314, 0.00415), 2)
  par <- avar[pt$par_free[pt$dest == "Lambda"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_lambda")

  # 5.3 Theta
  comp <- c(
    0.01820, 0.02698, 0.03034,
    0.01421, 0.02492, 0.03872
  )
  par <- avar[pt$par_free[pt$dest == "Theta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(0.00617, 0.00527)
  par <- avar[pt$par_free[pt$dest == "mu_eta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_eta")

  # 6.2 Variances
  comp <- c(0.03743, 0.05461)
  par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_var")

  # 7. Latent-manifest Covariances
  comp <- c(0.00928, 0.01120)
  par <- avar[pt$par_free[pt$dest == "Sigma_z_lv"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_cov")
})


# ---------------------------------------------------
# TEST 8 - one latent, two manifest covariates - two groups
# ---------------------------------------------------
test_that("two-group one latent, two manifest covariates Poisson", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ eta1 + z12 + z21",
    lv = list(eta1 = c("z41", "z42", "z43")),
    group = "treat",
    data = example01,
    family = "poisson",
    se = TRUE
  )
  # Converged?
  conv <- fit@fit@fit$convergence
  expect_equal(conv, 0)

  # Correct parameter estimates?
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 48)
  expect_equal(length(avar), 42)

  # LOG-LIKELIHOOD
  comp <- 11.58595
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(
    2.36960, -0.10848, 0.08538,
    2.57967, -0.10664, 0.06168
  )
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 2.2 gamma
  comp <- c(-0.08794, -0.02418)
  par <- pt$par[pt$dest == "gamma"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_gamma")

  # 3. Overdispersion parameter
  comp <- c(0, 0)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(1.35900, 3.91136, 1.39321, 4.00164)
  par <- pt$par[pt$dest == "mu_z"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

  # 4.2 Variances
  comp <- c(1.59255, 0.95087, 1.49459, 0.83742)
  par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

  # 4.3 Covariance
  comp <- c(-0.63469, -0.57679)
  par <- pt$par[pt$dest == "Sigma_z" & pt$type == "cov"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_cov")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0, -0.09006, -0.44041), 2)
  par <- pt$par[pt$dest == "nu"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(1, 1.27810, 1.35125), 2)
  par <- pt$par[pt$dest == "Lambda"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_lambda")

  # 5.3 Theta
  comp <- c(
    1.52390, 1.46958, 1.44107,
    1.36535, 1.53183, 1.06577
  )
  par <- pt$par[pt$dest == "Theta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(1.58494, 1.64547)
  par <- pt$par[pt$dest == "mu_eta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_eta")

  # 6.2 Variances
  comp <- c(1.86770, 2.14256)
  par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_var")

  # 7. Latent-Manifest Covariances
  comp <- c(0.50107, -0.26665, 0.69155, -0.38644)
  par <- pt$par[pt$dest == "Sigma_z_lv"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_eta_cov")


  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(
    0.00798, 0.00023, 0.00035,
    0.00690, 0.00018, 0.00030
  )
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 2.2 gamma
  comp <- c(0.00023, 0.00013)
  par <- avar[pt$par_free[pt$dest == "gamma"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_gamma")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(0.00371, 0.00222, 0.00323, 0.00185)
  par <- avar[pt$par_free[pt$dest == "mu_z"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

  # 4.2 Variances
  comp <- c(0.01189, 0.00423, 0.01039, 0.00321)
  par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")

  # 4.3 Covarianes
  comp <- c(0.00449, 0.00368)
  par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "cov"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0.01280, 0.01500), 2)
  par <- avar[pt$par_free[pt$dest == "nu"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(0.00310, 0.00390), 2)
  par <- avar[pt$par_free[pt$dest == "Lambda"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_lambda")

  # 5.3 Theta
  comp <- c(
    0.01783, 0.02565, 0.02852,
    0.01414, 0.02406, 0.03561
  )
  par <- avar[pt$par_free[pt$dest == "Theta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(0.00625, 0.00537)
  par <- avar[pt$par_free[pt$dest == "mu_eta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_eta")

  # 6.2 Variances
  comp <- c(0.03870, 0.05667)
  par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_var")

  # 7. Latent-manifest Covariances
  comp <- c(0.00952, 0.00536, 0.01146, 0.00591)
  par <- avar[pt$par_free[pt$dest == "Sigma_z_lv"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_cov")
})


test_that("two-group one latent, two manifest covariate negative binomial", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ eta1 + z12 + z21",
    lv = list(eta1 = c("z41", "z42", "z43")),
    group = "treat",
    data = example01,
    family = "nbinom",
    se = TRUE
  )

  # Converged?
  conv <- fit@fit@fit$convergence
  expect_equal(conv, 0)

  # Correct parameter estimates?
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 48)
  expect_equal(length(avar), 44)

  # LOG-LIKELIHOOD
  comp <- 11.49911
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(
    2.34399, -0.11310, 0.08741,
    2.58951, -0.10876, 0.05917
  )
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 2.2 gamma
  comp <- c(-0.07169, -0.02223)
  par <- pt$par[pt$dest == "gamma"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_gamma")

  # 3. Overdispersion parameter
  comp <- c(14.52583, 19.11102)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(1.35808, 3.91184, 1.39380, 4.00131)
  par <- pt$par[pt$dest == "mu_z"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

  # 4.2 Variances
  comp <- c(1.59218, 0.95074, 1.49559, 0.83776)
  par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

  # 4.3 Covariance
  comp <- c(-0.63448, -0.57736)
  par <- pt$par[pt$dest == "Sigma_z" & pt$type == "cov"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_cov")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0, -0.10279, -0.46896), 2)
  par <- pt$par[pt$dest == "nu"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(1, 1.28584, 1.36858), 2)
  par <- pt$par[pt$dest == "Lambda"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_lambda")

  # 5.3 Theta
  comp <- c(
    1.55152, 1.45768, 1.36277,
    1.38046, 1.53710, 1.03101
  )
  par <- pt$par[pt$dest == "Theta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(1.58175, 1.64726)
  par <- pt$par[pt$dest == "mu_eta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_eta")

  # 6.2 Variances
  comp <- c(1.84751, 2.12428)
  par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_var")

  # 7. Latent-Manifest Covariances
  comp <- c(0.49468, -0.26216, 0.68968, -0.38534)
  par <- pt$par[pt$dest == "Sigma_z_lv"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_eta_cov")


  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(
    0.01312, 0.00037, 0.00058,
    0.01155, 0.00030, 0.00050
  )
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 2.2 gamma
  comp <- c(0.00030, 0.00020)
  par <- avar[pt$par_free[pt$dest == "gamma"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_gamma")

  # 3. Overdispersion parameter
  comp <- c(5.62609, 10.32605)
  par <- avar[pt$par_free[pt$dest == "overdis"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(0.00370, 0.00222, 0.00322, 0.00184)
  par <- avar[pt$par_free[pt$dest == "mu_z"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

  # 4.2 Variances
  comp <- c(0.01188, 0.00423, 0.01041, 0.00322)
  par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")

  # 4.3 Covariances
  comp <- c(0.00449, 0.00368)
  par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "cov"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")


  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0.01297, 0.01581), 2)
  par <- avar[pt$par_free[pt$dest == "nu"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(0.00314, 0.00418), 2)
  par <- avar[pt$par_free[pt$dest == "Lambda"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_lambda")

  # 5.3 Theta
  comp <- c(
    0.01820, 0.02711, 0.03058,
    0.01426, 0.02493, 0.03876
  )
  par <- avar[pt$par_free[pt$dest == "Theta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(0.00617, 0.00527)
  par <- avar[pt$par_free[pt$dest == "mu_eta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_eta")

  # 6.2 Variances
  comp <- c(0.03746, 0.05544)
  par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_var")

  # 7. Latent-manifest Covariances
  comp <- c(0.00927, 0.00525, 0.01135, 0.00583)
  par <- avar[pt$par_free[pt$dest == "Sigma_z_lv"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_cov")
})


# ---------------------------------------------------
# TEST 9 - two latent, one manifest covariates - two groups
# ---------------------------------------------------
test_that("two-group two latent, one manifest covariates Poisson", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ eta1 + eta2 + z11",
    lv = list(
      eta1 = c("z21", "z22"),
      eta2 = c("z41", "z42", "z43")
    ),
    group = "treat",
    data = example01,
    family = "poisson",
    se = TRUE
  )
  # Converged?
  conv <- fit@fit@fit$convergence
  expect_equal(conv, 0)

  # Correct parameter estimates?
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 70)
  expect_equal(length(avar), 50)

  # LOG-LIKELIHOOD
  comp <- 12.55314
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09358)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(0.82196, -0.02026, 0.02869, 0.06399)
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 2.2 gamma
  comp <- c(0.42366, -0.03738, 0.61331, 0.03109)
  par <- pt$par[pt$dest == "gamma"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_gamma")

  # 3. Overdispersion parameter
  comp <- c(0, 0)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(1.59222, 1.55037)
  par <- pt$par[pt$dest == "mu_z"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

  # 4.2 Variances
  comp <- c(1.68941, 1.37974)
  par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0, 1.45202, 0, -0.10394, -0.45037), 2)
  par <- pt$par[pt$dest == "nu"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(1, 0.72596, rep(0, 5), 1, 1.28661, 1.35739), 2)
  par <- pt$par[pt$dest == "Lambda"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_lambda")

  # 5.3 Theta
  comp <- c(
    0.50781, 0.33968, 1.54610, 1.41844, 1.41326,
    0.52843, 0.37247, 1.36521, 1.53639, 1.05942
  )
  par <- pt$par[pt$dest == "Theta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(3.93821, 1.57983, 3.97780, 1.64804)
  par <- pt$par[pt$dest == "mu_eta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_eta")

  # 6.2 Variances
  comp <- c(0.45290, 1.88845, 0.30698, 2.14891)
  par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_var")

  # 6.3 Covariances
  comp <- c(-0.34781, -0.41573)
  par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "cov"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_cov")

  # 7. Latent-Manifest Covariances
  comp <- c(-0.50655, 0.47312, -0.45868, 0.64211)
  par <- pt$par[pt$dest == "Sigma_z_lv"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_eta_cov")


  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(0.13553, 0.00061, 0.31977, 0.00150)
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 2.2 gamma
  comp <- c(0.00629, 0.00042, 0.01502, 0.00050)
  par <- avar[pt$par_free[pt$dest == "gamma"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_gamma")

  # # 3. Overdispersion parameter
  # comp <- c(8.39126, 6.93554)
  # par <- avar[pt$par_free[pt$dest == "overdis"]] |> round(5)
  # expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(0.00392, 0.00294)
  par <- avar[pt$par_free[pt$dest == "mu_z"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

  # 4.2 Variances
  comp <- c(0.01346, 0.00892)
  par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")


  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0.04251, 0.01284, 0.01435), 2)
  par <- avar[pt$par_free[pt$dest == "nu"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(0.00268, 0.00310, 0.00365), 2)
  par <- avar[pt$par_free[pt$dest == "Lambda"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_lambda")

  # 5.3 Theta
  comp <- c(
    0.00349, 0.00117, 0.01784, 0.02584, 0.02760,
    0.00239, 0.00098, 0.01372, 0.02304, 0.02940
  )
  par <- avar[pt$par_free[pt$dest == "Theta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(0.00198, 0.00617, 0.00156, 0.00523)
  par <- avar[pt$par_free[pt$dest == "mu_eta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_eta")

  # 6.2 Variances
  comp <- c(0.00487, 0.04076, 0.00267, 0.05716)
  par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_var")

  # 6.3 Covariances
  comp <- c(0.00486, 0.00500)
  par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "cov"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_cov")

  # 7. Latent-manifest Covariances
  comp <- c(0.00384, 0.01017, 0.00270, 0.01060)
  par <- avar[pt$par_free[pt$dest == "Sigma_z_lv"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_cov")
})


test_that("two-group two latent, one manifest covariate negative binomial", {
  skip_on_cran()
  fit <- countreg(
    forml = "dv ~ eta1 + eta2 + z11",
    lv = list(
      eta1 = c("z21", "z22"),
      eta2 = c("z41", "z42", "z43")
    ),
    group = "treat",
    data = example01,
    family = "negbin",
    se = TRUE
  )
  # Converged?
  conv <- fit@fit@fit$convergence
  expect_equal(conv, 0)

  # Correct parameter estimates?
  pt <- fit@partable
  avar <- diag(fit@vcov)
  expect_equal(length(pt$par), 70)
  expect_equal(length(avar), 52)

  # LOG-LIKELIHOOD
  comp <- 12.50705
  par <- fit@fit@fit$objective
  expect_equal(par, comp, tolerance = 1e-5, label = "logl")

  # PARAMETER
  # 1. Group weight
  comp <- c(6.05912, 6.09357)
  par <- pt$par[pt$dest == "groupw"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(1.86183, -0.08019, 1.92432, -0.04937)
  par <- pt$par[pt$dest == "beta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

  # 2.2 gamma
  comp <- c(0.19708, -0.05938, 0.20441, -0.01444)
  par <- pt$par[pt$dest == "gamma"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_gamma")

  # 3. Overdispersion parameter
  comp <- c(16.07668, 19.65088)
  par <- pt$par[pt$dest == "overdis"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(1.59228, 1.55054)
  par <- pt$par[pt$dest == "mu_z"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

  # 4.2 Variances
  comp <- c(1.68916, 1.38302)
  par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0, 1.50289, 0, -0.10129, -0.44922), 2)
  par <- pt$par[pt$dest == "nu"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(1, 0.71290, rep(0, 5), 1, 1.28496, 1.35665), 2)
  par <- pt$par[pt$dest == "Lambda"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_lambda")

  # 5.3 Theta
  comp <- c(
    0.40008, 0.29254, 1.54578, 1.41983, 1.40770,
    0.39539, 0.31152, 1.36214, 1.54148, 1.05999
  )
  par <- pt$par[pt$dest == "Theta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(3.93729, 1.58028, 3.98050, 1.64902)
  par <- pt$par[pt$dest == "mu_eta"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_eta")

  # 6.2 Variances
  comp <- c(0.56083, 1.89159, 0.44179, 2.18016)
  par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "var"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_var")

  # 6.3 Covariances
  comp <- c(-0.34963, -0.42280)
  par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "cov"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_cov")

  # 7. Latent-Manifest Covariances
  comp <- c(-0.51018, 0.47128, -0.46614, 0.65297)
  par <- pt$par[pt$dest == "Sigma_z_lv"]
  expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_eta_cov")


  # STANDARD ERRORS
  # 1. Group weight
  comp <- c(0.00234, 0.00226)
  par <- avar[pt$par_free[pt$dest == "groupw"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

  # 2. Regression coefficient
  # 2.1 beta
  comp <- c(0.04109, 0.00038, 0.05660, 0.00046)
  par <- avar[pt$par_free[pt$dest == "beta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

  # 2.2 gamma
  comp <- c(0.00190, 0.00033, 0.00260, 0.00024)
  par <- avar[pt$par_free[pt$dest == "gamma"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_gamma")

  # 3. Overdispersion parameter
  comp <- c(8.68960, 13.23318)
  par <- avar[pt$par_free[pt$dest == "overdis"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

  # 4. Manifest Covariate Parameters
  # 4.1 Means
  comp <- c(0.00392, 0.00293)
  par <- avar[pt$par_free[pt$dest == "mu_z"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

  # 4.2 Variances
  comp <- c(0.01345, 0.00897)
  par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")


  # 5. Measurement Model
  # 5.1 nu
  comp <- rep(c(0.04897, 0.01268, 0.01420), 2)
  par <- avar[pt$par_free[pt$dest == "nu"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_nu")

  # 5.2 Lambda
  comp <- rep(c(0.00310, 0.00305, 0.00360), 2)
  par <- avar[pt$par_free[pt$dest == "Lambda"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_lambda")

  # 5.3 Theta
  comp <- c(
    0.00293, 0.00097, 0.01771, 0.02552, 0.02740,
    0.00221, 0.00081, 0.01349, 0.02321, 0.02914
  )
  par <- avar[pt$par_free[pt$dest == "Theta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_theta")

  # 6. Latent Covariate Parameters
  # 6.1 Means
  comp <- c(0.00205, 0.00618, 0.00164, 0.00513)
  par <- avar[pt$par_free[pt$dest == "mu_eta"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_eta")

  # 6.2 Variances
  comp <- c(0.00497, 0.04050, 0.00355, 0.05295)
  par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "var"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_var")

  # 6.3 Covariances
  comp <- c(0.00518, 0.00536)
  par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "cov"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_cov")

  # 7. Latent-manifest Covariances
  comp <- c(0.00409, 0.01013, 0.00301, 0.01050)
  par <- avar[pt$par_free[pt$dest == "Sigma_z_lv"]]
  expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_cov")
})


# ---------------------------------------------------
# TEST 10 - one latent variable - one group
# ---------------------------------------------------
test_that("one latent variable in one group - Poisson", {
  skip("Not finished")
  fit <- countreg(
    forml = "dv ~ eta",
    lv = list(eta = c("z41", "z42", "z43")),
    group = NULL,
    data = example01,
    family = "poisson",
    se = TRUE
  )
  par <- fit@partable$par
  comp <- c(
    6.769642, 2.72414, -0.109569, 0, 1, -0.054137,
    1.256452, -0.343201, 1.29337, 1.61712, 2.004875, 0,
    1.393845, 1.523669, 1.430018
  )
  expect_equal(length(par), 15)
  expect_equal(par, comp, tolerance = 1e-5)
})
