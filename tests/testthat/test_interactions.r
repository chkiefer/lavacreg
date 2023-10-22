# ---------------------------------------------------
# TEST 1 - two interacting manifest variables - one group
# ---------------------------------------------------
test_that("two interacting manifest variables in one group - Poisson", {
    skip("Not finished")
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

    # STANDARD ERRORS
    # 1. Group weight
    comp <- c(0.00234, 0.00226)
    par <- avar[pt$par_free[pt$dest == "groupw"]]
    expect_equal(par, comp, tolerance = 1e-5, label = "se_groupw")
})
