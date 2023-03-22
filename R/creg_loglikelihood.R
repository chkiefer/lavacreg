#' Fit the lavacreg model
#'
#' A wrapper for starting values, optimizing loglik and
#' computation of standard errors
#'
#'  @param object A lavacreg object
#'
#' @importFrom stats nlminb
#' @importFrom pracma hessian
#' @keywords internal
#' @noRd
creg_fit_model <- function(object) {
  input <- object@input
  family <- input@family
  no_groups <- input@no_groups
  no_lv <- input@no_lv
  silent <- input@silent
  se <- input@se

  dataobj <- object@dataobj
  datalist <- dataobj@datalist

  pt <- object@partable

  # Computation of start values
  if (no_lv > 0L) {
    # Start of Timer
    if (!silent) {
      cat("Computing starting values...")
      time_start <- Sys.time()
    }

    # Actual computation of start values if latent variables are involved
    # TODO: verschieben nach countreg
    x_start <- creg_starts_lv(object)

    # End of Timer
    if (!silent) {
      time_diff <- Sys.time() - time_start
      units(time_diff) <- "secs"
      cat("done. Took:", round(time_diff, 1), "s\n")
    }
  } else {
    # Use default starting values
    x_start <- matrix(pt$par[pt$par_free > 0L], nrow = 1)
  }


  # Constraints currently only used for measurement invariance
  # in latent variables
  # TODO: das kann schon viel frueher passieren und in eigenem Schritt
  if (no_lv > 0L & no_groups >= 2L) {
    # Get matrices for constraints
    tmp <- creg_constraints(pt)
    dataobj@eq_constraints_Q2 <- tmp$Q2
    dataobj@con_jac <- tmp$con_jac

    # Transform x-vector to "shorter" version
    x_start <- x_start %*% dataobj@eq_constraints_Q2
  }

  ##################
  # ESTIMATION PART
  ##################
  # Objective function
  objective_function <- function(x) {
    # CHECK for missing values in x
    if (anyNA(x)) {
      return(+Inf)
    }

    # Specify x as column-vector (for orientation in matrix computations)
    x <- matrix(x, ncol = 1)

    # If equality constraints exist, re-expand x to full length
    if (no_lv > 0L & no_groups >= 2L) {
      x <- as.numeric(dataobj@eq_constraints_Q2 %*% x)
    }

    # Save currents iteration values to partable (connect to meaning)
    pt$par[pt$par_free > 0L] <- x

    # NO variances equal or smaller to 0
    if (any(pt$par[pt$dest == "lv_grid" & pt$type == "var"] <= 0)) {
      return(+Inf)
    }

    # Create modellist (pre-compute all kinds of stuff for model estimation)
    # TODO: skip modellist and hand pt directly over to C++ likelihood function
    modellist <- creg_modellist(
      pt = pt,
      dataobj = dataobj,
      family = family,
      input = input
    )

    # Call loglikelihood function to compute actual likelihood
    obj <- creg_loglikelihood_function(datalist, modellist)
    return(obj)
  }

  # Start of Timer
  if (!silent) {
    cat("Fitting the model...")
    time_start <- Sys.time()
  }

  # Pass start values x and objective function to optimizer
  # TODO: maybe allow for different optimizers
  fit <- nlminb(x_start, objective_function,
    control = list(
      rel.tol = 1e-6,
      eval.max = 500,
      iter.max = 300
    )
  )

  # End of Timer
  if (!silent) {
    time_diff <- Sys.time() - time_start
    units(time_diff) <- "secs"
    cat("done. Took:", round(time_diff, 1), "s\n")
  }

  #####################################
  # Here starts the standard error part
  ####################################
  if (!fit$convergence & se) {
    if (!silent) {
      cat("Computing standard errors...")
      time_start <- Sys.time()
    }
    information <- hessian(objective_function, fit$par)
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
      solve(information) / sum(object@input@n_cell),
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

  if (no_lv & no_groups >= 2) {
    pt$par[pt$par_free > 0L] <- as.numeric(
      dataobj@eq_constraints_Q2 %*% fit$par
    )
    if (!is.null(vcov_fit)) {
      vcov_fit <-
        dataobj@eq_constraints_Q2 %*% vcov_fit %*% t(dataobj@eq_constraints_Q2)
    }
  } else {
    pt$par[pt$par_free > 0L] <- fit$par
  }

  object@fit <- list(
    fit = fit,
    vcov_fit = vcov_fit,
    pt = pt
  )
  return(object)
}


#' lavacreg log-likelihood function
#'
#' Computes the likelihood function using the C++ function for
#' group-conditional likelihoods
#'
#' @param datalist Datalist
#' @param modellist Modellist
#'
#' @importFrom stats dpois
#' @keywords internal
#' @noRd
creg_loglikelihood_function <- function(datalist, modellist) {
  kappas <- modellist$groupw
  n_cell <- modellist$n_cell
  no_groups <- length(kappas)
  family <- modellist$family

  obj_outgroup <- sum(dpois(n_cell, exp(kappas), log = TRUE))

  obj_ingroups <- mapply(function(data, modellist_g) {
    muy <- modellist_g$muy
    sigmayw <- modellist_g$sigmayw
    muwz <- modellist_g$muwz
    sigmaz <- modellist_g$sigmaz
    ghweight <- modellist_g$ghweight
    detvarz <- modellist_g$detvarz
    dims <- modellist_g$dims

    if (any(!is.na(sigmaz))) {
      if (any(diag(solve(sigmaz)) <= 0)) {
        return(-Inf)
      }
    }
    if (any(sigmayw[-1] <= 0)) {
      return(-Inf)
    }
    if (family == "nbinom" & sigmayw[1] <= 0) {
      return(-Inf)
    }

    obj_i <- compute_groupcond_logl(
      x = data, muy = muy, sigmayw = sigmayw, muwz = muwz,
      sigmaz = sigmaz, ghweight = ghweight, detvarz = detvarz,
      dims = dims
    )
    if (is.na(obj_i)) {
      return(-Inf)
    }
    return(obj_i)
  }, data = datalist, modellist_g = modellist$modellist_g, SIMPLIFY = TRUE)

  obj <- -(obj_outgroup + sum(obj_ingroups)) / sum(n_cell)
  if (is.na(obj)) {
    return(+Inf)
  }
  return(obj)
}
