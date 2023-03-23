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
creg_model_estimate <- function(object) {
  input <- object@input
  family <- input@family
  no_groups <- input@no_groups
  no_lv <- input@no_lv
  silent <- input@silent
  se <- input@se

  datalist <- object@datalist
  constraints <- object@constraints
  gh_grid <- object@gh_grid
  pt <- object@partable
  x_start <- object@x_start


  # Objective function
  objective_function <- function(x) {
    # CHECK for missing values in x
    if (anyNA(x)) {
      return(+Inf)
    }

    # Specify x as column-vector (for orientation in matrix computations)
    x <- matrix(x, ncol = 1)

    # If equality constraints exist, re-expand x to full length
    if (constraints@con_logical) {
      x <- as.numeric(constraints@eq_constraints_Q2 %*% x)
    }

    # Save current iteration values to partable (connect to meaning)
    pt$par[pt$par_free > 0L] <- x

    # NO variances equal or smaller to 0
    if (any(pt$par[pt$dest == "lv_grid" & pt$type == "var"] <= 0)) {
      return(+Inf)
    }

    # Create modellist (pre-compute all kinds of stuff for model estimation)
    # TODO: skip modellist and hand pt directly over to C++ likelihood function
    modellist <- creg_modellist(
      pt = pt,
      datalist = datalist,
      gh_grid = gh_grid,
      family = family,
      input = input
    )

    # Call loglikelihood function to compute actual likelihood
    obj <- creg_model_objective(datalist, modellist)
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

  # Save parameters back to partable
  if (constraints@con_logical) {
    pt$par[pt$par_free > 0L] <- as.numeric(
      constraints@eq_constraints_Q2 %*% fit$par
    )
  } else {
    pt$par[pt$par_free > 0L] <- fit$par
  }

  # End of Timer
  if (!silent) {
    time_diff <- Sys.time() - time_start
    units(time_diff) <- "secs"
    cat("done. Took:", round(time_diff, 1), "s\n")
  }

  # Save fit and partable back to object and return object
  object@fit <- new("creg_fit",
    fit = fit,
    objective = objective_function
  )
  object@partable <- pt

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
creg_model_objective <- function(datalist, modellist) {
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
