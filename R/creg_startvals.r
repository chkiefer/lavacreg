#' Starting values for parameter estimation
#'
#' Derive starting values
#'
#'  @param object a lavacreg object
#'
#' @noRd
creg_startvals <- function(object) {
    no_lv <- object@input@no_lv
    silent <- object@input@silent
    pt <- object@partable
    constraints <- object@constraints

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

    # Constraints currently only used for measurement invariance
    # in latent variables
    if (constraints@con_logical) {
        # Transform x-vector to "shorter" version
        x_start <- x_start %*% constraints@eq_constraints_Q2
    }

    # End of Timer
    if (!silent) {
        time_diff <- Sys.time() - time_start
        units(time_diff) <- "secs"
        cat("done. Took:", round(time_diff, 1), "s\n")
    }

    return(x_start)
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
#' @noRd
creg_starts_lv <- function(object) {
    # TODO: integrate starting values for measurement model as well
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
    # covariances
    pt_starts <- fit_starts@partable

    # STARTING VALUES FOR
    # 1. Group weight
    pt$par[pt$dest == "groupw"] <- pt_starts$par[pt_starts$dest == "groupw"]

    # 2. Regression coefficients
    pt$par[pt$dest == "regcoef"] <- pt_starts$par[pt_starts$dest == "regcoef"]

    # 3. Covariate Means
    pt$par[!is.na(pt$type) & pt$type == "mean"] <-
        pt_starts$par[!is.na(pt_starts$type) & pt_starts$type == "mean"]

    # 4. Covariate Variances
    pt$par[!is.na(pt$type) & pt$type == "var"] <-
        pt_starts$par[!is.na(pt_starts$type) & pt_starts$type == "var"]

    # 5. Covariate Covariances
    pt$par[!is.na(pt$type) & (pt$type == "cov" | pt$type == "cov_z_lv")] <-
        pt_starts$par[!is.na(pt_starts$type) & pt_starts$type == "cov"]

    # Return starting values
    x_start_lv <- pt$par[pt$par_free > 0L] |> matrix(nrow = 1)
    return(x_start_lv)
}
