# ---------------------------------------------------
# TEST 1 - two interacting manifest variables - one group
# ---------------------------------------------------
test_that("two interacting manifest variables in one group - Poisson", {
    # skip("Not finished")
    fit <- countreg(
        forml = "dv ~ z12*z21",
        group = NULL,
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
    expect_equal(length(pt$par), 11)
    expect_equal(length(avar), 10)

    # LOG-LIKELIHOOD
    comp <- 5.989974
    par <- fit@fit@fit$objective
    expect_equal(par, comp, tolerance = 1e-5, label = "logl")

    # PARAMETER
    # 1. Group weight
    comp <- c(6.76964)
    par <- pt$par[pt$dest == "groupw"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(2.30403, -0.09048, 0.09862)
    par <- pt$par[pt$dest == "beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

    # 2.2 Beta
    comp <- c(-0.006605)
    par <- pt$par[pt$dest == "Beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_Beta")

    # 3. Overdispersion parameter
    comp <- c(0)
    par <- pt$par[pt$dest == "overdis"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

    # 4. Manifest Covariate Parameters
    # 4.1 Means
    comp <- c(1.38462, 3.95274)
    par <- pt$par[pt$dest == "mu_z"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

    # 4.2 Variances
    comp <- c(1.53481, 0.89243)
    par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

    # 4.3 Covariances
    comp <- c(-0.59969)
    par <- pt$par[pt$dest == "Sigma_z" & pt$type == "cov"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_cov")

    # STANDARD ERRORS
    # 1. Group weight
    comp <- c(0.001148)
    par <- avar[pt$par_free[pt$dest == "groupw"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(0.00700, 0.00121, 0.00035)
    par <- avar[pt$par_free[pt$dest == "beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

    # 2.2 Beta
    comp <- c(0.00007)
    par <- avar[pt$par_free[pt$dest == "Beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_Beta")

    # # 3. Overdispersion parameter
    # comp <- c(8.39126, 6.93554)
    # par <- avar[pt$par_free[pt$dest == "overdis"]] |> round(5)
    # expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

    # 4. Manifest Covariate Parameters
    # 4.1 Means
    comp <- c(0.00176, 0.00102)
    par <- avar[pt$par_free[pt$dest == "mu_z"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

    # 4.2 Variances
    comp <- c(0.00541, 0.00183)
    par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")

    # 4.3 Covariances
    comp <- c(0.00199)
    par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "cov"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")
})


test_that(
    "two interacting manifest variables in one group - negative binomial",
    {
        # skip("Not finished")
        fit <- countreg(
            forml = "dv ~ z12*z21",
            group = NULL,
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
        expect_equal(length(pt$par), 11)
        expect_equal(length(avar), 11)

        # LOG-LIKELIHOOD
        comp <- 5.86041
        par <- fit@fit@fit$objective
        expect_equal(par, comp, tolerance = 1e-5, label = "logl")

        # PARAMETER
        # 1. Group weight
        comp <- c(6.76964)
        par <- pt$par[pt$dest == "groupw"]
        expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

        # 2. Regression coefficient
        # 2.1 beta
        comp <- c(2.32888, -0.10141, 0.09320)
        par <- pt$par[pt$dest == "beta"]
        expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

        # 2.2 Beta
        comp <- c(-0.00422)
        par <- pt$par[pt$dest == "Beta"]
        expect_equal(par, comp, tolerance = 1e-5, label = "par_Beta")

        # 3. Overdispersion parameter
        comp <- c(13.2323)
        par <- pt$par[pt$dest == "overdis"]
        expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

        # 4. Manifest Covariate Parameters
        # 4.1 Means
        comp <- c(1.38462, 3.95274)
        par <- pt$par[pt$dest == "mu_z"]
        expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

        # 4.2 Variances
        comp <- c(1.53481, 0.89243)
        par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
        expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

        # 4.3 Covariances
        comp <- c(-0.59969)
        par <- pt$par[pt$dest == "Sigma_z" & pt$type == "cov"]
        expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_cov")

        # STANDARD ERRORS
        # 1. Group weight
        comp <- c(0.001148)
        par <- avar[pt$par_free[pt$dest == "groupw"]]
        expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

        # 2. Regression coefficient
        # 2.1 beta
        comp <- c(0.01351, 0.00221, 0.00068)
        par <- avar[pt$par_free[pt$dest == "beta"]]
        expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

        # 2.2 Beta
        comp <- c(0.00013)
        par <- avar[pt$par_free[pt$dest == "Beta"]]
        expect_equal(par, comp, tolerance = 1e-5, label = "se_Beta")

        # 3. Overdispersion parameter
        comp <- c(1.85621)
        par <- avar[pt$par_free[pt$dest == "overdis"]]
        expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

        # 4. Manifest Covariate Parameters
        # 4.1 Means
        comp <- c(0.00176, 0.00102)
        par <- avar[pt$par_free[pt$dest == "mu_z"]]
        expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

        # 4.2 Variances
        comp <- c(0.00541, 0.00183)
        par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
        expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")

        # 4.3 Covariances
        comp <- c(0.00199)
        par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "cov"]]
        expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")
    }
)


# ---------------------------------------------------
# TEST 2 - two interacting manifest variables - two group
# ---------------------------------------------------
test_that("two interacting manifest variables in two groups - Poisson", {
    # skip("Not finished")
    fit <- countreg(
        forml = "dv ~ z12*z21",
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
    expect_equal(length(pt$par), 22)
    expect_equal(length(avar), 20)

    # LOG-LIKELIHOOD
    comp <- 5.91916
    par <- fit@fit@fit$objective
    expect_equal(par, comp, tolerance = 1e-5, label = "logl")

    # PARAMETER
    # 1. Group weight
    comp <- c(6.05913, 6.09357)
    par <- pt$par[pt$dest == "groupw"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(
        2.26994, -0.15421, 0.08528,
        2.36591, -0.02328, 0.10459
    )
    par <- pt$par[pt$dest == "beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

    # 2.2 Beta
    comp <- c(0.00589, -0.02279)
    par <- pt$par[pt$dest == "Beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_Beta")

    # 3. Overdispersion parameter
    comp <- c(0, 0)
    par <- pt$par[pt$dest == "overdis"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

    # 4. Manifest Covariate Parameters
    # 4.1 Means
    comp <- c(1.36294, 3.90926, 1.40556, 3.99474)
    par <- pt$par[pt$dest == "mu_z"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

    # 4.2 Variances
    comp <- c(1.59049, 0.95031, 1.48013, 0.83293)
    par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

    # 4.3 Covariances
    comp <- c(-0.63361, -0.56873)
    par <- pt$par[pt$dest == "Sigma_z" & pt$type == "cov"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_cov")

    # STANDARD ERRORS
    # 1. Group weight
    comp <- c(0.00234, 0.00226)
    par <- avar[pt$par_free[pt$dest == "groupw"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(0.01356, 0.00238, 0.00068, 0.01497, 0.00265, 0.00074)
    par <- avar[pt$par_free[pt$dest == "beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

    # 2.2 Beta
    comp <- c(0.00014, 0.00015)
    par <- avar[pt$par_free[pt$dest == "Beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_Beta")

    # # 3. Overdispersion parameter
    # comp <- c(8.39126, 6.93554)
    # par <- avar[pt$par_free[pt$dest == "overdis"]] |> round(5)
    # expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

    # 4. Manifest Covariate Parameters
    # 4.1 Means
    comp <- c(0.00372, 0.00222, 0.00334, 0.00188)
    par <- avar[pt$par_free[pt$dest == "mu_z"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

    # 4.2 Variances
    comp <- c(0.01182, 0.00422, 0.00989, 0.00313)
    par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")

    # 4.3 Covariances
    comp <- c(0.00447, 0.00351)
    par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "cov"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")
})


test_that("two interacting manifest variables in two groups - NB", {
    # skip("Not finished")
    fit <- countreg(
        forml = "dv ~ z12*z21",
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
    expect_equal(length(pt$par), 22)
    expect_equal(length(avar), 22)

    # LOG-LIKELIHOOD
    comp <- 5.82045
    par <- fit@fit@fit$objective
    expect_equal(par, comp, tolerance = 1e-5, label = "logl")

    # PARAMETER
    # 1. Group weight
    comp <- c(6.05913, 6.09357)
    par <- pt$par[pt$dest == "groupw"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(
        2.29131, -0.16483,  0.08037,
        2.37998, -0.02839,  0.10173
    )
    par <- pt$par[pt$dest == "beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

    # 2.2 Beta
    comp <- c(0.00849, -0.02191)
    par <- pt$par[pt$dest == "Beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_Beta")

    # 3. Overdispersion parameter
    comp <- c(12.99607, 19.05656)
    par <- pt$par[pt$dest == "overdis"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

    # 4. Manifest Covariate Parameters
    # 4.1 Means
    comp <- c(1.36293, 3.90927, 1.40557, 3.99473)
    par <- pt$par[pt$dest == "mu_z"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

    # 4.2 Variances
    comp <- c(1.59050, 0.95030, 1.48011, 0.83293)
    par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

    # 4.3 Covariances
    comp <- c(-0.63361, -0.56872)
    par <- pt$par[pt$dest == "Sigma_z" & pt$type == "cov"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_cov")

    # STANDARD ERRORS
    # 1. Group weight
    comp <- c(0.00234, 0.00226)
    par <- avar[pt$par_free[pt$dest == "groupw"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(0.02471, 0.00398, 0.00125, 0.02605, 0.00446, 0.00129)
    par <- avar[pt$par_free[pt$dest == "beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

    # 2.2 Beta
    comp <- c(0.00024, 0.00026)
    par <- avar[pt$par_free[pt$dest == "Beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_Beta")

    # 3. Overdispersion parameter
    comp <- c(3.88956, 10.26889)
    par <- avar[pt$par_free[pt$dest == "overdis"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

    # 4. Manifest Covariate Parameters
    # 4.1 Means
    comp <- c(0.00372, 0.00222, 0.00334, 0.00188)
    par <- avar[pt$par_free[pt$dest == "mu_z"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

    # 4.2 Variances
    comp <- c(0.01182, 0.00422, 0.00989, 0.00313)
    par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")

    # 4.3 Covariances
    comp <- c(0.00447, 0.00351)
    par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "cov"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")
})


# ---------------------------------------------------
# TEST 3 - two interacting and additional manifest variables - one group
# ---------------------------------------------------
test_that("two interacting and additional manifest variables in one group - Poisson", {
    # skip("Not finished")
    fit <- countreg(
        forml = "dv ~ z12*z21 + z11",
        group = NULL,
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
    expect_equal(length(pt$par), 18)
    expect_equal(length(avar), 15)

    # LOG-LIKELIHOOD
    comp <- 7.345426
    par <- fit@fit@fit$objective
    expect_equal(par, comp, tolerance = 1e-5, label = "logl")

    # PARAMETER
    # 1. Group weight
    comp <- c(6.76964)
    par <- pt$par[pt$dest == "groupw"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(2.40624, -0.05666, 0.08613, -0.06369)
    par <- pt$par[pt$dest == "beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

    # 2.2 Beta
    comp <- c(-0.00680, 0.00000, 0.00000)
    par <- pt$par[pt$dest == "Beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_Beta")

    # 3. Overdispersion parameter
    comp <- c(0)
    par <- pt$par[pt$dest == "overdis"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

    # 4. Manifest Covariate Parameters
    # 4.1 Means
    comp <- c(1.38462, 3.95273, 1.57980)
    par <- pt$par[pt$dest == "mu_z"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

    # 4.2 Variances
    comp <- c(1.53482, 0.89243, 1.52101)
    par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

    # 4.3 Covariances
    comp <- c(-0.59969, 0.94137, -0.48630)
    par <- pt$par[pt$dest == "Sigma_z" & pt$type == "cov"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_cov")

    # STANDARD ERRORS
    # 1. Group weight
    comp <- c(0.001148)
    par <- avar[pt$par_free[pt$dest == "groupw"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(0.00723, 0.00124, 0.00035, 0.00010)
    par <- avar[pt$par_free[pt$dest == "beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

    # 2.2 Beta
    comp <- c(0.00007)
    par <- avar[pt$par_free[pt$dest == "Beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_Beta")

    # # 3. Overdispersion parameter
    # comp <- c(8.39126, 6.93554)
    # par <- avar[pt$par_free[pt$dest == "overdis"]] |> round(5)
    # expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

    # 4. Manifest Covariate Parameters
    # 4.1 Means
    comp <- c(0.00176, 0.00102, 0.00175)
    par <- avar[pt$par_free[pt$dest == "mu_z"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

    # 4.2 Variances
    comp <- c(0.00541, 0.00183, 0.00531)
    par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")

    # 4.3 Covariances
    comp <- c(0.00199, 0.003698, 0.001830)
    par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "cov"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")
})

test_that("two interacting and additional manifest variables in one group - NB", {
    # skip("Not finished")
    fit <- countreg(
        forml = "dv ~ z12*z21 + z11",
        group = NULL,
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
    expect_equal(length(pt$par), 18)
    expect_equal(length(avar), 16)

    # LOG-LIKELIHOOD
    comp <- 7.22667
    par <- fit@fit@fit$objective
    expect_equal(par, comp, tolerance = 1e-5, label = "logl")

    # PARAMETER
    # 1. Group weight
    comp <- c(6.76964)
    par <- pt$par[pt$dest == "groupw"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(2.42867, -0.06671, 0.08132, -0.06452)
    par <- pt$par[pt$dest == "beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

    # 2.2 Beta
    comp <- c(-0.00440, 0.00000, 0.00000)
    par <- pt$par[pt$dest == "Beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_Beta")

    # 3. Overdispersion parameter
    comp <- c(13.96401)
    par <- pt$par[pt$dest == "overdis"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

    # 4. Manifest Covariate Parameters
    # 4.1 Means
    comp <- c(1.38462, 3.95274, 1.57979)
    par <- pt$par[pt$dest == "mu_z"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

    # 4.2 Variances
    comp <- c(1.53481, 0.89243, 1.52103)
    par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

    # 4.3 Covariances
    comp <- c(-0.59969, 0.94137, -0.48631)
    par <- pt$par[pt$dest == "Sigma_z" & pt$type == "cov"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_cov")

    # STANDARD ERRORS
    # 1. Group weight
    comp <- c(0.001148)
    par <- avar[pt$par_free[pt$dest == "groupw"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(0.01365, 0.00222, 0.00067, 0.00019)
    par <- avar[pt$par_free[pt$dest == "beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

    # 2.2 Beta
    comp <- c(0.00013)
    par <- avar[pt$par_free[pt$dest == "Beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_Beta")

    # 3. Overdispersion parameter
    comp <- c(2.18779)
    par <- avar[pt$par_free[pt$dest == "overdis"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

    # 4. Manifest Covariate Parameters
    # 4.1 Means
    comp <- c(0.00176, 0.00102, 0.00175)
    par <- avar[pt$par_free[pt$dest == "mu_z"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

    # 4.2 Variances
    comp <- c(0.00541, 0.00183, 0.00531)
    par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")

    # 4.3 Covariances
    comp <- c(0.00199, 0.003698, 0.001830)
    par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "cov"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")
})


# ---------------------------------------------------
# TEST 4 - two interacting manifest and additional latent variables - one group
# ---------------------------------------------------
test_that("two interacting manifest and additional latent variables in one group - Poisson", {
    # skip("Not finished")
    fit <- countreg(
        forml = "dv ~ z12*z21 + eta",
        lv = list(eta = c("z41", "z42", "z43")),
        group = NULL,
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
    expect_equal(length(pt$par), 25)
    expect_equal(length(avar), 22)

    # LOG-LIKELIHOOD
    comp <- 11.66712
    par <- fit@fit@fit$objective
    expect_equal(par, comp, tolerance = 1e-5, label = "logl")

    # PARAMETER
    # 1. Group weight
    comp <- c(6.76964)
    par <- pt$par[pt$dest == "groupw"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(2.38657, -0.07559, 0.09206)
    par <- pt$par[pt$dest == "beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

    # 2.2 gamma
    comp <- c(-0.04848)
    par <- pt$par[pt$dest == "gamma"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_gamma")

    # 2.3 Beta
    comp <- c(-0.00651)
    par <- pt$par[pt$dest == "Beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_Beta")

    # 3. Overdispersion parameter
    comp <- c(0)
    par <- pt$par[pt$dest == "overdis"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

    # 4. Manifest Covariate Parameters
    # 4.1 Means
    comp <- c(1.37621, 3.95727)
    par <- pt$par[pt$dest == "mu_z"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_z")

    # 4.2 Variances
    comp <- c(1.53908, 0.89368) # mismatch 1
    par <- pt$par[pt$dest == "Sigma_z" & pt$type == "var"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_var")

    # 4.3 Covariances
    comp <- c(-0.60200) # mismatch
    par <- pt$par[pt$dest == "Sigma_z" & pt$type == "cov"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_cov")

    # 5. Measurement Model
    # 5.1 nu
    comp <- c(0, -0.09690, -0.41836) # mismatch 2
    par <- pt$par[pt$dest == "nu"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_nu")

    # 5.2 Lambda
    comp <- c(1, 1.28191, 1.33711)
    par <- pt$par[pt$dest == "Lambda"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_lambda")

    # 5.3 Theta
    comp <- c(1.43791, 1.47194, 1.29449)
    par <- pt$par[pt$dest == "Theta" & pt$type == "var"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_theta")

    # 6. Latent Covariate Parameters
    # 6.1 Means
    comp <- c(1.61464)
    par <- pt$par[pt$dest == "mu_eta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_eta")

    # 6.2 Variances
    comp <- c(1.97684)
    par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "var"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_var")

    # 7. Latent-Manifest Covariances
    comp <- c(0.58781, -0.31690)
    par <- pt$par[pt$dest == "Sigma_z_lv"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_z_eta_cov")

    # STANDARD ERRORS
    # 1. Group weight
    comp <- c(0.001148)
    par <- avar[pt$par_free[pt$dest == "groupw"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(0.00727, 0.00122, 0.00035)
    par <- avar[pt$par_free[pt$dest == "beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

    # 2.2 gamma
    comp <- c(0.00008)
    par <- avar[pt$par_free[pt$dest == "gamma"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_gamma")

    # 2.2 Beta
    comp <- c(0.00007)
    par <- avar[pt$par_free[pt$dest == "Beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_Beta")

    # # 3. Overdispersion parameter
    # comp <- c(8.39126, 6.93554)
    # par <- avar[pt$par_free[pt$dest == "overdis"]] |> round(5)
    # expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

    # 4. Manifest Covariate Parameters
    # 4.1 Means
    comp <- c(0.00174, 0.00102)
    par <- avar[pt$par_free[pt$dest == "mu_z"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_z")

    # 4.2 Variances
    comp <- c(0.00548, 0.00184)
    par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "var"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")

    # 4.3 Covariances
    comp <- c(0.00201)
    par <- avar[pt$par_free[pt$dest == "Sigma_z" & pt$type == "cov"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_z_var")

    # 5. Measurement Model
    # 5.1 nu
    comp <- c(0.01294, 0.01428)
    par <- avar[pt$par_free[pt$dest == "nu"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_nu")

    # 5.2 Lambda
    comp <- c(0.00316, 0.00365)
    par <- avar[pt$par_free[pt$dest == "Lambda"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_lambda")

    # 5.3 Theta
    comp <- c(0.00852, 0.01391, 0.01714)
    par <- avar[pt$par_free[pt$dest == "Theta" & pt$type == "var"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_theta")

    # 6. Latent Covariate Parameters
    # 6.1 Means
    comp <- c(0.00363)
    par <- avar[pt$par_free[pt$dest == "mu_eta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_eta")

    # 6.2 Variances
    comp <- c(0.02873)
    par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "var"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_var")

    # 7. Latent-manifest Covariances
    comp <- c(0.00513, 0.00272)
    par <- avar[pt$par_free[pt$dest == "Sigma_z_lv"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_cov")
})



# ---------------------------------------------------
# TEST 5 - two interacting latent variables - one group
# ---------------------------------------------------
test_that("two interacting manifest and additional latent variables in one group - Poisson", {
    # skip("Not finished")
    fit <- countreg(
        forml = "dv ~ eta1*eta2",
        lv = list(
            eta1 = c("z11", "z12"),
            eta2 = c("z41", "z42", "z43")
        ),
        group = NULL,
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
    expect_equal(length(pt$par), 31)
    expect_equal(length(avar), 21)

    # LOG-LIKELIHOOD
    comp <- 11.77966
    par <- fit@fit@fit$objective
    expect_equal(par, comp, tolerance = 1e-5, label = "logl")

    # PARAMETER
    # 1. Group weight
    comp <- c(6.76964)
    par <- pt$par[pt$dest == "groupw"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(3.00449)
    par <- pt$par[pt$dest == "beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

    # 2.2 gamma
    comp <- c(-0.31990, 0.04171)
    par <- pt$par[pt$dest == "gamma"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_gamma")

    # 2.4 Gamma
    comp <- c(-0.01770)
    par <- pt$par[pt$dest == "Gamma"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_Beta")

    # 3. Overdispersion parameter
    comp <- c(0)
    par <- pt$par[pt$dest == "overdis"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

    # 5. Measurement Model
    # 5.1 nu
    comp <- c(0, -0.30762, 0, -0.11155, -0.40829)
    par <- pt$par[pt$dest == "nu"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_nu")

    # 5.2 Lambda
    comp <- c(1.07118, 1.29081, 1.33096)
    par <- pt$par[pt$dest == "Lambda" & pt$par_free]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_lambda")

    # 5.3 Theta
    comp <- c(0.74783, 0.66189, 1.45928, 1.42279, 1.32088)
    par <- pt$par[pt$dest == "Theta" & pt$type == "var"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_theta")

    # 6. Latent Covariate Parameters
    # 6.1 Means
    comp <- c(1.57099, 1.61952)
    par <- pt$par[pt$dest == "mu_eta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_eta")

    # 6.2 Variances
    comp <- c(0.76438, 1.99461)
    par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "var"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_var")

    # 6.3 Covariances
    comp <- c(0.54075)
    par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "cov"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_cov")

    # STANDARD ERRORS
    # 1. Group weight
    comp <- c(0.001148)
    par <- avar[pt$par_free[pt$dest == "groupw"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(0.00185)
    par <- avar[pt$par_free[pt$dest == "beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

    # 2.2 gamma
    comp <- c(0.00097, 0.00067)
    par <- avar[pt$par_free[pt$dest == "gamma"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_gamma")

    # 2.4 Gamma
    comp <- c(0.00015)
    par <- avar[pt$par_free[pt$dest == "Gamma"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_Beta")

    # # 3. Overdispersion parameter
    # comp <- c(8.39126, 6.93554)
    # par <- avar[pt$par_free[pt$dest == "overdis"]] |> round(5)
    # expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

    # 5. Measurement Model
    # 5.1 nu
    comp <- c(0.00996, 0.01270, 0.01304)
    par <- avar[pt$par_free[pt$dest == "nu"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_nu")

    # 5.2 Lambda
    comp <- c(0.00329, 0.00307, 0.00317)
    par <- avar[pt$par_free[pt$dest == "Lambda"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_lambda")

    # 5.3 Theta
    comp <- c(0.00268, 0.00277, 0.00809, 0.01377, 0.01440)
    par <- avar[pt$par_free[pt$dest == "Theta" & pt$type == "var"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_theta")

    # 6. Latent Covariate Parameters
    # 6.1 Means
    comp <- c(0.00171, 0.00377)
    par <- avar[pt$par_free[pt$dest == "mu_eta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_eta")

    # 6.2 Variances
    comp <- c(0.00522, 0.02907)
    par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "var"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_var")

    # 6.3 Covariances
    comp <- c(0.00409)
    par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "cov"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_cov")
})

# ---------------------------------------------------
# TEST 6 - one latent and one manifest interaction - one group
# ---------------------------------------------------
test_that("two interacting manifest and additional latent variables in one group - Poisson", {
    skip("Not finished")
    fit <- countreg(
        forml = "dv ~ eta1*eta2 + z21*z22",
        lv = list(
            eta1 = c("z11", "z12"),
            eta2 = c("z41", "z42", "z43")
        ),
        group = NULL,
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
    expect_equal(length(pt$par), 31)
    expect_equal(length(avar), 21)

    # LOG-LIKELIHOOD
    comp <- 11.77966
    par <- fit@fit@fit$objective
    expect_equal(par, comp, tolerance = 1e-5, label = "logl")

    # PARAMETER
    # 1. Group weight
    comp <- c(6.76964)
    par <- pt$par[pt$dest == "groupw"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(3.00449)
    par <- pt$par[pt$dest == "beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

    # 2.2 gamma
    comp <- c(-0.31990, 0.04171)
    par <- pt$par[pt$dest == "gamma"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_gamma")

    # 2.4 Gamma
    comp <- c(-0.01770)
    par <- pt$par[pt$dest == "Gamma"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_Beta")

    # 3. Overdispersion parameter
    comp <- c(0)
    par <- pt$par[pt$dest == "overdis"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

    # 5. Measurement Model
    # 5.1 nu
    comp <- c(0, -0.30762, 0, -0.11155, -0.40829)
    par <- pt$par[pt$dest == "nu"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_nu")

    # 5.2 Lambda
    comp <- c(1.07118, 1.29081, 1.33096)
    par <- pt$par[pt$dest == "Lambda" & pt$par_free]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_lambda")

    # 5.3 Theta
    comp <- c(0.74783, 0.66189, 1.45928, 1.42279, 1.32088)
    par <- pt$par[pt$dest == "Theta" & pt$type == "var"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_theta")

    # 6. Latent Covariate Parameters
    # 6.1 Means
    comp <- c(1.57099, 1.61952)
    par <- pt$par[pt$dest == "mu_eta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_eta")

    # 6.2 Variances
    comp <- c(0.76438, 1.99461)
    par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "var"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_var")

    # 6.3 Covariances
    comp <- c(0.54075)
    par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "cov"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_cov")

    # STANDARD ERRORS
    # 1. Group weight
    comp <- c(0.001148)
    par <- avar[pt$par_free[pt$dest == "groupw"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(0.00185)
    par <- avar[pt$par_free[pt$dest == "beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

    # 2.2 gamma
    comp <- c(0.00097, 0.00067)
    par <- avar[pt$par_free[pt$dest == "gamma"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_gamma")

    # 2.4 Gamma
    comp <- c(0.00015)
    par <- avar[pt$par_free[pt$dest == "Gamma"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_Beta")

    # # 3. Overdispersion parameter
    # comp <- c(8.39126, 6.93554)
    # par <- avar[pt$par_free[pt$dest == "overdis"]] |> round(5)
    # expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

    # 5. Measurement Model
    # 5.1 nu
    comp <- c(0.00996, 0.01270, 0.01304)
    par <- avar[pt$par_free[pt$dest == "nu"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_nu")

    # 5.2 Lambda
    comp <- c(0.00329, 0.00307, 0.00317)
    par <- avar[pt$par_free[pt$dest == "Lambda"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_lambda")

    # 5.3 Theta
    comp <- c(0.00268, 0.00277, 0.00809, 0.01377, 0.01440)
    par <- avar[pt$par_free[pt$dest == "Theta" & pt$type == "var"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_theta")

    # 6. Latent Covariate Parameters
    # 6.1 Means
    comp <- c(0.00171, 0.00377)
    par <- avar[pt$par_free[pt$dest == "mu_eta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_eta")

    # 6.2 Variances
    comp <- c(0.00522, 0.02907)
    par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "var"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_var")

    # 6.3 Covariances
    comp <- c(0.00409)
    par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "cov"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_cov")
})

# ---------------------------------------------------
# TEST 6 - one latent and one manifest interaction - one group
# ---------------------------------------------------
test_that("mixed interaction - Poisson", {
    skip("Not finished")
    fit <- countreg(
        forml = "dv ~ eta1*z22",
        lv = list(eta1 = c("z11", "z12")),
        group = NULL,
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
    expect_equal(length(pt$par), 31)
    expect_equal(length(avar), 21)

    # LOG-LIKELIHOOD
    comp <- 11.77966
    par <- fit@fit@fit$objective
    expect_equal(par, comp, tolerance = 1e-5, label = "logl")

    # PARAMETER
    # 1. Group weight
    comp <- c(6.76964)
    par <- pt$par[pt$dest == "groupw"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(3.00449)
    par <- pt$par[pt$dest == "beta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_beta")

    # 2.2 gamma
    comp <- c(-0.31990, 0.04171)
    par <- pt$par[pt$dest == "gamma"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_gamma")

    # 2.4 Gamma
    comp <- c(-0.01770)
    par <- pt$par[pt$dest == "Gamma"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_Beta")

    # 3. Overdispersion parameter
    comp <- c(0)
    par <- pt$par[pt$dest == "overdis"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_overdis")

    # 5. Measurement Model
    # 5.1 nu
    comp <- c(0, -0.30762, 0, -0.11155, -0.40829)
    par <- pt$par[pt$dest == "nu"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_nu")

    # 5.2 Lambda
    comp <- c(1.07118, 1.29081, 1.33096)
    par <- pt$par[pt$dest == "Lambda" & pt$par_free]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_lambda")

    # 5.3 Theta
    comp <- c(0.74783, 0.66189, 1.45928, 1.42279, 1.32088)
    par <- pt$par[pt$dest == "Theta" & pt$type == "var"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mm_theta")

    # 6. Latent Covariate Parameters
    # 6.1 Means
    comp <- c(1.57099, 1.61952)
    par <- pt$par[pt$dest == "mu_eta"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_mu_eta")

    # 6.2 Variances
    comp <- c(0.76438, 1.99461)
    par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "var"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_var")

    # 6.3 Covariances
    comp <- c(0.54075)
    par <- pt$par[pt$dest == "Sigma_eta" & pt$type == "cov"]
    expect_equal(par, comp, tolerance = 1e-5, label = "par_sig_eta_cov")

    # STANDARD ERRORS
    # 1. Group weight
    comp <- c(0.001148)
    par <- avar[pt$par_free[pt$dest == "groupw"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")

    # 2. Regression coefficient
    # 2.1 beta
    comp <- c(0.00185)
    par <- avar[pt$par_free[pt$dest == "beta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_beta")

    # 2.2 gamma
    comp <- c(0.00097, 0.00067)
    par <- avar[pt$par_free[pt$dest == "gamma"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_gamma")

    # 2.4 Gamma
    comp <- c(0.00015)
    par <- avar[pt$par_free[pt$dest == "Gamma"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_Beta")

    # # 3. Overdispersion parameter
    # comp <- c(8.39126, 6.93554)
    # par <- avar[pt$par_free[pt$dest == "overdis"]] |> round(5)
    # expect_equal(par, comp, tolerance = 1e-5, label = "se_overdis")

    # 5. Measurement Model
    # 5.1 nu
    comp <- c(0.00996, 0.01270, 0.01304)
    par <- avar[pt$par_free[pt$dest == "nu"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_nu")

    # 5.2 Lambda
    comp <- c(0.00329, 0.00307, 0.00317)
    par <- avar[pt$par_free[pt$dest == "Lambda"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_lambda")

    # 5.3 Theta
    comp <- c(0.00268, 0.00277, 0.00809, 0.01377, 0.01440)
    par <- avar[pt$par_free[pt$dest == "Theta" & pt$type == "var"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mm_theta")

    # 6. Latent Covariate Parameters
    # 6.1 Means
    comp <- c(0.00171, 0.00377)
    par <- avar[pt$par_free[pt$dest == "mu_eta"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_mu_eta")

    # 6.2 Variances
    comp <- c(0.00522, 0.02907)
    par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "var"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_var")

    # 6.3 Covariances
    comp <- c(0.00409)
    par <- avar[pt$par_free[pt$dest == "Sigma_eta" & pt$type == "cov"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_sig_eta_cov")
})
