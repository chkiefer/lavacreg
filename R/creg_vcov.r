creg_vcov <- function(object, information_only = FALSE) {
    fit <- object@fit@fit
    objective_function <- object@fit@objective
    se <- object@input@se
    silent <- object@input@silent
    n_cell <- object@input@n_cell
    constraints <- object@constraints
    pt <- object@partable
    par_est <- pt$par[pt$par_free > 0]

    if (constraints@con_logical) {
        # Transform x-vector to "shorter" version
        par_est <- par_est %*% constraints@eq_constraints_Q2
    }

    gh_grid <- object@gh_grid
    datalist <- object@datalist
    family <- object@input@family
    input <- object@input

    #####################################
    # Here starts the standard error part
    ####################################
    if (!fit$convergence & se) {
        if (!silent) {
            cat("Computing standard errors...")
            time_start <- Sys.time()
        }
        # browser()
        information <- hessian(objective_function, par_est)

        if (information_only) {
            if (constraints@con_logical) {
                if (!is.null(information)) {
                    information <-
                        constraints@eq_constraints_Q2 %*%
                        information %*%
                        t(constraints@eq_constraints_Q2)
                }
            }
            if (!silent) {
                time_diff <- Sys.time() - time_start
                units(time_diff) <- "secs"
                cat("done. Took:", round(time_diff, 1), "s\n")
            }
            return(information)
        }

        eigvals <- eigen(information,
            symmetric = TRUE,
            only.values = TRUE
        )$values
        if (any(eigvals < -1 * .Machine$double.eps^(3 / 4))) {
            warning(
                "lavacreg WARNING: information matrix is not positive definite;
        the model may not be identified"
            )
        }
        vcov_fit <- try(
            solve(information) / sum(n_cell),
            silent = TRUE
        )
        if (!silent) {
            time_diff <- Sys.time() - time_start
            units(time_diff) <- "secs"
            cat("done. Took:", round(time_diff, 1), "s\n")
        }
    } else if (fit$convergence & se) {
        vcov_fit <- NULL
        warning(
            "lavacreg warning: Estimation did not converge.
      Standard errors are not computed."
        )
    } else {
        (
            vcov_fit <- NULL
        )
    }

    if (constraints@con_logical) {
        if (!is.null(vcov_fit)) {
            vcov_fit <-
                constraints@eq_constraints_Q2 %*%
                vcov_fit %*%
                t(constraints@eq_constraints_Q2)
        }
    }
    return(vcov_fit)
}
