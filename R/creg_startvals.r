#' Starting values for parameter estimation
#'
#' Derive starting values
#'
#'  @param object a lavacreg object
#'
#' @noRd
creg_startvals <- function(object) {
    input <- object@input
    no_lv <- input@no_lv
    silent <- input@silent
    pt <- object@partable
    # constraints <- object@constraints

    # Start of Timer
    if (!silent) {
        cat("Computing starting values...")
        time_start <- Sys.time()
    }

    if (no_lv > 0L) {
        # Actual computation of start values if latent variables are involved
        x_start <- creg_starts_lv(object)
    } else {
        # Use default starting values
        x_start <- pt$par[pt$par_free > 0L] |> matrix(nrow = 1)
    }

    # Save starting values to partable
    pt$start <- pt$par
    pt$start[pt$par_free > 0L] <- x_start

    # End of Timer
    if (!silent) {
        time_diff <- Sys.time() - time_start
        units(time_diff) <- "secs"
        cat("done. Took:", round(time_diff, 1), "s\n")
    }

    return(pt)
}





#' Starting values for latent variable countreg models
#'
#' This function computes starting values for countreg models
#' involving latent covariates. Starting values are computed
#' as follows:
#' 1. Mean scores for latent variable indicators are computed
#' 2. A count regression with these mean scores are estimated
#' 3. Parameters of manifest model are returned as starting values
#'
#' This procedure is not perfect and is meant to reduce the time spent
#' in numerical integration, by getting as close to the plausible estimates
#' as possible without the MML estimation.
#'
#'  @param object a lavacreg object
#'
#' @importFrom stats var
#' @noRd
creg_starts_lv <- function(object) {
    # TODO: integrate starting values for measurement model as well
    input <- object@input
    dvname <- input@dvname
    lvnames <- input@lvnames
    ovnames <- input@ovnames
    cvnames <- input@cvnames
    intnames <- input@intnames
    no_int_z <- input@no_int_z
    no_int_lv <- input@no_int_lv
    no_int_z_lv <- input@no_int_z_lv
    no_lv <- input@no_lv
    n_cell <- input@n_cell
    lvlist <- input@lvlist
    groupname <- input@groupname
    data <- input@data
    cfa <- input@cfa
    pt <- object@partable


    # Group part
    pt$par[pt$dest == "groupw"] <- log(n_cell)


    ## MEASUREMENT MODEL PART
    # Adapted from lavaan
    COV <- data[ovnames] |> var()

    pt$par[pt$dest == "Theta"] <- COV |>
        diag() * 0.5

    for (lv in names(lvlist)) {
        ovs <- lvlist[[lv]]
        sub_COV <- COV[ovs, ovs]
        fabin <- lavaan:::lav_cfa_1fac_fabin(
            sub_COV,
            std.lv = FALSE,
            lambda.only = TRUE,
            method = "fabin3"
        )
        pt$par[pt$dest == "Lambda" & pt$lhs == lv & pt$rhs %in% ovs] <- fabin$lambda
    }

    ov_free_names <- pt$lhs[pt$dest == "nu" & pt$par_free]
    ov_free_means <- data[ov_free_names] |> colMeans()
    pt$par[pt$dest == "nu" & pt$par_free] <- ov_free_means


    # 1. Regression intercept
    # mu_Y <- unlist(data[dvname]) |> mean()
    # pt$par[pt$dest == "beta" & pt$lhs %in% dvname] <- log(mu_Y)


    if (!cfa) {
        # browser()
        # Compute (manifest) mean scores for latent variables
        # for (lv in names(lvlist)) {
        #     data[lv] <- rowMeans(data[lvlist[[lv]]], na.rm = TRUE)
        # }

        Lambda_B <- pt$par[pt$dest == "Lambda" & pt$group == 1] |> matrix(ncol = no_lv)
        Theta_B_inv <- pt$par[pt$dest == "Theta" & pt$group == 1] |>
            diag() |>
            solve()

        A_B <- solve(t(Lambda_B) %*% Theta_B_inv %*% Lambda_B) %*% t(Lambda_B) %*% Theta_B_inv
        fac_scores <- as.matrix(data[ovnames]) %*% t(A_B)
        colnames(fac_scores) <- names(lvlist)
        fac_scores <- as.data.frame(fac_scores)
        for (lv in names(lvlist)) {
            data[lv] <- fac_scores[lv]
        }

        # Dataset of all model variables including mean score variables
        d_ov_starts <- data[c(dvname, groupname, cvnames, lvnames)]

        # Model formula for "all-observed" case
        forml <- object@input@forml

        # alter Code war wegen Reihenfolge, oder?
        # forml <- paste(
        #     dvname,
        #     "~",
        #     paste(cvnames, collapse = "+"),
        #     "+",
        #     paste(lvnames, collapse = "+")
        # )

        # Run the "all-observed" model through countreg
        fit_starts <- countreg(
            forml = forml,
            lv = NULL,
            group = groupname,
            data = d_ov_starts,
            family = object@input@family,
            silent = TRUE,
            se = FALSE
        )

        # Extract "manifest" results for coefficients, means, variance, and
        # covariances
        pt_starts <- fit_starts@partable
        # browser()

        # STARTING VALUES FOR
        # 1. Group weight
        pt$par[pt$dest == "groupw"] <- pt_starts$par[pt_starts$dest == "groupw"]

        # 2. Regression coefficients
        pt$par[pt$dest == "beta"] <-
            pt_starts$par[
                pt_starts$dest == "beta" & pt_starts$rhs %in% c(1, cvnames)
            ]
        pt$par[pt$dest == "gamma"] <-
            pt_starts$par[pt_starts$dest == "beta" & pt_starts$rhs %in% lvnames]


        # Interactions
        if (no_int_z) {
            model_z_int <- apply(intnames$z, 1, paste, collapse = ":")
            model_z_int2 <- apply(
                matrix(intnames$z[, c(2, 1)], ncol = 2), 1, paste,
                collapse = ":"
            )

            pt$par[pt$dest == "Beta" & pt$par_free] <-
                pt_starts$par[
                    pt_starts$dest == "Beta" &
                        pt_starts$par_free &
                        pt_starts$rhs %in% c(model_z_int, model_z_int2)
                ]
        }


        if (no_int_lv) {
            model_lv_int <- apply(intnames$lv, 1, paste, collapse = ":")
            model_lv_int2 <- apply(
                matrix(intnames$lv[, c(2, 1)], ncol = 2), 1, paste,
                collapse = ":"
            )

            pt$par[pt$dest == "Gamma" & pt$par_free] <-
                pt_starts$par[
                    pt_starts$dest == "Beta" &
                        pt_starts$par_free &
                        pt_starts$rhs %in% c(model_lv_int, model_lv_int2)
                ]
        }

        if (no_int_z_lv) {
            model_z_lv_int <- apply(intnames$z_lv, 1, paste, collapse = ":")
            model_z_lv_int2 <- apply(
                matrix(intnames$z_lv[, c(2, 1)], ncol = 2), 1, paste,
                collapse = ":"
            )

            pt$par[pt$dest == "Omega" & pt$par_free] <-
                pt_starts$par[
                    pt_starts$dest == "Beta" &
                        pt_starts$par_free &
                        pt_starts$rhs %in% c(model_z_lv_int, model_z_lv_int2)
                ]
        }






        # 3. Covariate Means
        pt$par[pt$dest == "mu_z"] <-
            pt_starts$par[pt_starts$dest == "mu_z" & pt_starts$lhs %in% cvnames]
        pt$par[pt$dest == "mu_eta"] <-
            pt_starts$par[pt_starts$dest == "mu_z" & pt_starts$lhs %in% lvnames]

        # 4. Covariate Variances
        pt$par[(pt$dest == "Sigma_z") & pt$type == "var"] <-
            pt_starts$par[
                pt_starts$dest == "Sigma_z" & pt_starts$type == "var" & pt_starts$lhs %in% cvnames
            ]
        pt$par[(pt$dest == "Sigma_eta") & pt$type == "var"] <-
            pt_starts$par[
                pt_starts$dest == "Sigma_z" & pt_starts$type == "var" & pt_starts$lhs %in% lvnames
            ]

        # 5. Covariate Covariances
        pt$par[(pt$dest == "Sigma_z") & pt$type == "cov"] <-
            pt_starts$par[
                pt_starts$dest == "Sigma_z" & pt_starts$type == "cov" & pt_starts$lhs %in% cvnames & pt_starts$rhs %in% cvnames
            ]
        pt$par[(pt$dest == "Sigma_eta") & pt$type == "cov"] <-
            pt_starts$par[
                pt_starts$dest == "Sigma_z" & pt_starts$type == "cov" & pt_starts$lhs %in% lvnames & pt_starts$rhs %in% lvnames
            ]

        pt$par[(pt$dest == "Sigma_z_lv")] <-
            pt_starts$par[
                pt_starts$dest == "Sigma_z" & pt_starts$type == "cov" & ((pt_starts$lhs %in% lvnames & pt_starts$rhs %in% cvnames) |
                    (pt_starts$lhs %in% cvnames & pt_starts$rhs %in% lvnames))
            ]
    }


    # Return starting values
    x_start_lv <- pt$par[pt$par_free > 0L] |> matrix(nrow = 1)
    return(x_start_lv)
}



creg_starts_to_partable <- function(x, constraints) {
    browser()
    x <- matrix(x, ncol = 1)

    # If equality constraints exist, re-expand x to full length
    if (constraints@con_logical) {
        x <- as.numeric(constraints@eq_constraints_Q2 %*% x)
    }
    x
}
