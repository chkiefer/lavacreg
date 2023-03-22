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

    # End of Timer
    if (!silent) {
        time_diff <- Sys.time() - time_start
        units(time_diff) <- "secs"
        cat("done. Took:", round(time_diff, 1), "s\n")
    }

    return(x_start)
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

    x_start_lv <- pt$par[pt$par_free > 0L] |> matrix(nrow = 1)

    return(x_start_lv)
}
