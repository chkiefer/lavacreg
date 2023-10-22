#' @title Create An Initial Parameter Table for Count Regression Model
#'
#' @description This function turns the input into a parameter table.
#' The paramater table (partable) serves three purposes:
#' 1. Translating the input into the required parameters
#' 2. The initial \code{par} column contains starting values for each free
#' parameter
#' 3. Later, it holds the values of each iteration within the fitting
#' process and connects them to model
#'
#' @details Details
#'
#' @param object An object of class \code{data.frame}, which contains an
#' overview of parameters to be estimated
#'
#' @importFrom utils combn
#' @keywords internal
#' @noRd
creg_partable <- function(input) {
    # Import required information from input object
    dvname <- input@dvname
    lvnames <- input@lvnames
    ovnames <- input@ovnames
    cvnames <- input@cvnames
    intnames <- input@intnames
    groupname <- input@groupname
    no_groups <- input@no_groups
    no_lv <- input@no_lv
    no_w <- input@no_w
    no_z <- input@no_z
    no_int_z <- input@no_int_z
    no_int_lv <- input@no_int_lv
    no_int_z_lv <- input@no_int_z_lv
    lv <- input@lvlist
    family <- input@family

    # Compute number of different covariances types
    no_cov <- no_lv + no_z
    no_z_lv_covariance <- no_z * no_lv
    no_z_covariance <- no_z * (no_z - 1) / 2
    no_lv_covariance <- no_lv * (no_lv - 1) / 2

    if (!length(groupname)) groupname <- ""

    covnames <- c(1, cvnames, lvnames)

    # Parts of the partable
    lhs <- NULL
    op <- NULL
    rhs <- NULL
    dest <- NULL
    type <- NULL
    group <- NULL
    par_free <- NULL
    par <- NULL

    # IDs and helpers
    par_free_id <- 0L

    # create partable
    for (i in 1:no_groups) {
        # Group weights
        lhs <- c(lhs, groupname)
        op <- c(op, "%")
        rhs <- c(rhs, "w")
        dest <- c(dest, "groupw")
        type <- c(type, NA)
        group <- c(group, i)
        par_free_id <- par_free_id + 1L
        par_free <- c(par_free, par_free_id)
        par <- c(par, 0)

        # Regression coefficients for Z
        for (j in 0:no_z) {
            lhs <- c(lhs, dvname)
            op <- c(op, "~")
            rhs <- c(rhs, covnames[j + 1])
            dest <- c(dest, "beta")
            type <- c(type, NA)
            group <- c(group, i)
            par_free_id <- par_free_id + 1L
            par_free <- c(par_free, par_free_id)
            par <- c(par, 0.0)
        }

        # Regression coefficients for eta
        if (no_lv) {
            for (j in 1:no_lv) {
                lhs <- c(lhs, dvname)
                op <- c(op, "~")
                rhs <- c(rhs, covnames[j + 1 + no_z])
                dest <- c(dest, "gamma")
                type <- c(type, NA)
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0.0)
            }
        }

        # Regression coefficients for interactions of Zs
        if (no_int_z) {
            names_z_int <- combn(cvnames, 2) |> apply(2, paste, collapse = ":")
            model_z_int <- apply(intnames$z, 1, paste, collapse = ":")
            model_z_int2 <- apply(
                matrix(intnames$z[, c(2, 1)], ncol = 2), 1, paste,
                collapse = ":"
            )
            for (j in names_z_int) {
                lhs <- c(lhs, dvname)
                op <- c(op, "~")
                rhs <- c(rhs, j)
                dest <- c(dest, "Beta")
                type <- c(type, NA)
                group <- c(group, i)
                if (j %in% model_z_int | j %in% model_z_int2) {
                    par_free_id <- par_free_id + 1L
                    par_free <- c(par_free, par_free_id)
                    par <- c(par, 0.0)
                } else {
                    par_free <- c(par_free, 0)
                    par <- c(par, 0.0)
                }
            }
        }

        # Regression coefficients for interactions of LVs
        if (no_int_lv) {
            names_lv_int <- combn(lvnames, 2) |> apply(2, paste, collapse = ":")
            model_lv_int <- apply(intnames$lv, 1, paste, collapse = ":")
            model_lv_int2 <- apply(
                matrix(intnames$lv[, c(2, 1)], ncol = 2), 1, paste,
                collapse = ":"
            )
            for (j in names_lv_int) {
                lhs <- c(lhs, dvname)
                op <- c(op, "~")
                rhs <- c(rhs, j)
                dest <- c(dest, "Gamma")
                type <- c(type, NA)
                group <- c(group, i)
                if (j %in% model_lv_int | j %in% model_lv_int2) {
                    par_free_id <- par_free_id + 1L
                    par_free <- c(par_free, par_free_id)
                    par <- c(par, 0.0)
                } else {
                    par_free <- c(par_free, 0)
                    par <- c(par, 0.0)
                }
            }
        }

        # Regression coefficients for interactions of Zs and LVs
        if (no_int_z_lv) {
            names_z_lv_int <- expand.grid(cvnames, lvnames) |> apply(2, paste, collapse = ":")
            model_z_lv_int <- apply(intnames$z_lv, 1, paste, collapse = ":")
            model_z_lv_int2 <- apply(
                matrix(intnames$z_lv[, c(2, 1)], ncol = 2), 1, paste,
                collapse = ":"
            )
            for (j in names_z_lv_int) {
                lhs <- c(lhs, dvname)
                op <- c(op, "~")
                rhs <- c(rhs, j)
                dest <- c(dest, "Omega")
                type <- c(type, NA)
                group <- c(group, i)
                if (j %in% model_z_lv_int | j %in% model_z_lv_int2) {
                    par_free_id <- par_free_id + 1L
                    par_free <- c(par_free, par_free_id)
                    par <- c(par, 0.0)
                } else {
                    par_free <- c(par_free, 0)
                    par <- c(par, 0.0)
                }
            }
        }

        # Overdispersion
        lhs <- c(lhs, dvname)
        op <- c(op, "~~")
        rhs <- c(rhs, dvname)
        dest <- c(dest, "overdis")
        type <- c(type, NA)
        group <- c(group, i)
        if (family == "poisson") {
            par_free <- c(par_free, 0)
            par <- c(par, 0)
        } else {
            par_free_id <- par_free_id + 1L
            par_free <- c(par_free, par_free_id)
            par <- c(par, 1)
        }




        # Measurement model
        no_ind_before <- 0L
        no_ind_after <- no_w

        if (no_lv) {
            for (j in 1:no_lv) {
                no_ind <- length(lv[[j]])
                no_ind_after <- no_ind_after - no_ind

                if (no_ind_before) {
                    for (l in 1:no_ind_before) {
                        lhs <- c(lhs, lvnames[j])
                        op <- c(op, "=~")
                        rhs <- c(rhs, ovnames[l])
                        dest <- c(dest, "Lambda")
                        type <- c(type, NA)
                        group <- c(group, i)
                        par_free <- c(par_free, 0)
                        par <- c(par, 0)
                    }
                }

                for (k in 1:no_ind) {
                    # Intercepts nu
                    lhs <- c(lhs, ovnames[no_ind_before + k])
                    op <- c(op, "~")
                    rhs <- c(rhs, 1)
                    dest <- c(dest, "nu")
                    type <- c(type, NA)
                    group <- c(group, i)

                    # Fix first intercept to zero
                    if (k != 1) {
                        par_free_id <- par_free_id + 1L
                        par_free <- c(par_free, par_free_id)
                        par <- c(par, 0)
                    } else {
                        par_free <- c(par_free, 0)
                        par <- c(par, 0)
                    }

                    # Factor loadings lambda
                    lhs <- c(lhs, lvnames[j])
                    op <- c(op, "=~")
                    rhs <- c(rhs, ovnames[no_ind_before + k])
                    dest <- c(dest, "Lambda")
                    type <- c(type, NA)
                    group <- c(group, i)

                    # Fix first loading to one
                    if (k != 1) {
                        par_free_id <- par_free_id + 1L
                        par_free <- c(par_free, par_free_id)
                        par <- c(par, 1)
                    } else {
                        par_free <- c(par_free, 0)
                        par <- c(par, 1)
                    }
                }

                if (no_ind_after) {
                    for (l in 1:no_ind_after) {
                        lhs <- c(lhs, lvnames[j])
                        op <- c(op, "=~")
                        rhs <- c(rhs, ovnames[no_ind_before + no_ind + l])
                        dest <- c(dest, "Lambda")
                        type <- c(type, NA)
                        group <- c(group, i)
                        par_free <- c(par_free, 0)
                        par <- c(par, 0)
                    }
                }

                no_ind_before <- no_ind_before + no_ind
            }

            # Measurement error variance
            for (j in 1:no_w) {
                lhs <- c(lhs, ovnames[j])
                op <- c(op, "~~")
                rhs <- c(rhs, ovnames[j])

                dest <- c(dest, "Theta")
                type <- c(type, "var")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 1)
            }
        }

        # Means and variances of latent variables
        if (no_lv) {
            for (j in 1:no_lv) {
                lhs <- c(lhs, lvnames[j])
                op <- c(op, "~")
                rhs <- c(rhs, 1)
                dest <- c(dest, "mu_eta")
                type <- c(type, NA)
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0)

                lhs <- c(lhs, lvnames[j])
                op <- c(op, "~~")
                rhs <- c(rhs, lvnames[j])
                dest <- c(dest, "Sigma_eta")
                type <- c(type, "var")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 1)
            }
        }

        if (no_lv_covariance) {
            names_lv_cov <- combn(lvnames, 2)
            for (j in 1:no_lv_covariance) {
                lhs <- c(lhs, names_lv_cov[1, j])
                op <- c(op, "~~")
                rhs <- c(rhs, names_lv_cov[2, j])
                dest <- c(dest, "Sigma_eta")
                type <- c(type, "cov")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0)
            }
        }


        # Means and variances of manifest covariates
        if (no_z) {
            for (j in 1:no_z) {
                lhs <- c(lhs, cvnames[j])
                op <- c(op, "~")
                rhs <- c(rhs, 1)
                dest <- c(dest, "mu_z")
                type <- c(type, NA)
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0)

                lhs <- c(lhs, cvnames[j])
                op <- c(op, "~~")
                rhs <- c(rhs, cvnames[j])
                dest <- c(dest, "Sigma_z")
                type <- c(type, "var")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 1)
            }
        }

        if (no_z_covariance) {
            names_z_cov <- combn(cvnames, 2)
            for (j in 1:no_z_covariance) {
                lhs <- c(lhs, names_z_cov[1, j])
                op <- c(op, "~~")
                rhs <- c(rhs, names_z_cov[2, j])
                dest <- c(dest, "Sigma_z")
                type <- c(type, "cov")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0)
            }
        }


        # Covariances between manifest and latent covariates
        if (no_z_lv_covariance) {
            names_z_lv_cov <- expand.grid(
                lvnames,
                cvnames,
                stringsAsFactors = FALSE
            )
            for (j in 1:no_z_lv_covariance) {
                lhs <- c(lhs, names_z_lv_cov[j, 1])
                op <- c(op, "~~")
                rhs <- c(rhs, names_z_lv_cov[j, 2])
                dest <- c(dest, "Sigma_z_lv")
                type <- c(type, NA)
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0)
            }
        }
    }
    pt <- data.frame(lhs, op, rhs, dest, type, group, par_free, par)
    return(pt)
}
