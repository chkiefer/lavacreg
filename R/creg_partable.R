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
creg_create_partable <- function(object) {
    # Import required information from input object
    input <- object@input
    dvname <- input@dvname
    lvnames <- input@lvnames
    ovnames <- input@ovnames
    cvnames <- input@cvnames
    groupname <- input@groupname
    no_groups <- input@no_groups
    no_lv <- input@no_lv
    no_w <- input@no_w
    no_z <- input@no_z
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

        # Regression coefficients
        for (j in 0:no_cov) {
            lhs <- c(lhs, dvname)
            op <- c(op, "~")
            rhs <- c(rhs, covnames[j + 1])
            dest <- c(dest, "regcoef")
            type <- c(type, NA)
            group <- c(group, i)
            par_free_id <- par_free_id + 1L
            par_free <- c(par_free, par_free_id)
            par <- c(par, 0.0)
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
                        dest <- c(dest, "mm")
                        type <- c(type, "lambda")
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
                    dest <- c(dest, "mm")
                    type <- c(type, "nu")
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
                    dest <- c(dest, "mm")
                    type <- c(type, "lambda")
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
                        dest <- c(dest, "mm")
                        type <- c(type, "lambda")
                        group <- c(group, i)
                        par_free <- c(par_free, 0)
                        par <- c(par, 0)
                    }
                }

                no_ind_before <- no_ind_before + no_ind
            }
        }

        # Means and variances of latent variables
        if (no_lv) {
            for (j in 1:no_lv) {
                lhs <- c(lhs, lvnames[j])
                op <- c(op, "~")
                rhs <- c(rhs, 1)
                dest <- c(dest, "lv_grid")
                type <- c(type, "mean")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0)

                lhs <- c(lhs, lvnames[j])
                op <- c(op, "~~")
                rhs <- c(rhs, lvnames[j])
                dest <- c(dest, "lv_grid")
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
                dest <- c(dest, "lv_grid")
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
                dest <- c(dest, "z")
                type <- c(type, "mean")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0)

                lhs <- c(lhs, cvnames[j])
                op <- c(op, "~~")
                rhs <- c(rhs, cvnames[j])
                dest <- c(dest, "z")
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
                dest <- c(dest, "z")
                type <- c(type, "cov")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0)
            }
        }


        # Covariances between manifest and latent covariates
        if (no_z_lv_covariance) {
            names_z_lv_cov <- expand.grid(lvnames, cvnames, stringsAsFactors = FALSE)
            for (j in 1:no_z_lv_covariance) {
                lhs <- c(lhs, names_z_lv_cov[j, 1])
                op <- c(op, "~~")
                rhs <- c(rhs, names_z_lv_cov[j, 2])
                dest <- c(dest, "z")
                type <- c(type, "cov_z_lv")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0)
            }
        }


        # Overdispersion and measurement error variance
        if (no_w) {
            for (j in 0:no_w) {
                if (j) {
                    lhs <- c(lhs, ovnames[j])
                    op <- c(op, "~~")
                    rhs <- c(rhs, ovnames[j])
                } else {
                    lhs <- c(lhs, dvname)
                    op <- c(op, "~~")
                    rhs <- c(rhs, dvname)
                }
                dest <- c(dest, "sigmaw")
                type <- c(type, ifelse(j, "veps", "size"))
                group <- c(group, i)
                if (j) {
                    par_free_id <- par_free_id + 1L
                    par_free <- c(par_free, par_free_id)
                } else if (family == "poisson") {
                    par_free <- c(par_free, 0)
                } else {
                    par_free_id <- par_free_id + 1L
                    par_free <- c(par_free, par_free_id)
                }

                par <- c(par, ifelse(j, 1, ifelse(family == "poisson", 0, 1)))
            }
        } else {
            lhs <- c(lhs, dvname)
            op <- c(op, "~~")
            rhs <- c(rhs, dvname)
            dest <- c(dest, "sigmaw")
            type <- c(type, "size")
            group <- c(group, i)
            if (family == "poisson") {
                par_free <- c(par_free, 0)
                par <- c(par, 0)
            } else {
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 1)
            }
        }
    }
    pt <- data.frame(lhs, op, rhs, dest, type, group, par_free, par)
    return(pt)
}

#' Start values
#'
#' Compute starting values for a latent variable model
#'
#'  @param object a lavacreg object
#'  @param pt a parameter table with initial starting values
#'
#' @noRd
creg_starts_lv <- function(object) {
    input <- object@input
    dvname <- input@dvname
    cvnames <- input@cvnames
    lvnames <- input@lvnames
    lvlist <- input@lvlist
    groupname <- input@groupname
    data <- input@data
    pt <- object@partable

    # Compute (manifest) mean scores for latent variables
    for (lv in names(lvlist)) {
        data[lv] <- rowMeans(data[lvlist[[lv]]], na.rm = TRUE)
    }

    # Dataset of all model variables including mean score variables
    d_ov_starts <- data[c(dvname, groupname, cvnames, lvnames)]

    # Model formula for "all-observed" case
    forml <- paste(
        dvname,
        "~",
        paste(cvnames, collapse = "+"),
        "+",
        paste(lvnames, collapse = "+")
    )

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
    # covariances and return as start values for latent model
    pt_starts <- fit_starts@fit$pt
    pt$par[pt$dest == "groupw"] <- pt_starts$par[pt_starts$dest == "groupw"]
    pt$par[pt$dest == "regcoef"] <- pt_starts$par[pt_starts$dest == "regcoef"]
    pt$par[!is.na(pt$type) & pt$type == "mean"] <-
        pt_starts$par[!is.na(pt_starts$type) & pt_starts$type == "mean"]
    pt$par[!is.na(pt$type) & pt$type == "var"] <-
        pt_starts$par[!is.na(pt_starts$type) & pt_starts$type == "var"]
    pt$par[!is.na(pt$type) & (pt$type == "cov" | pt$type == "cov_z_lv")] <-
        pt_starts$par[!is.na(pt_starts$type) & pt_starts$type == "cov"]

    x_start_lv <- pt$par[pt$par_free > 0L]
    return(x_start_lv)
}
