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
creg_model_estimate <- function(object, fit.model = TRUE) {
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
  x_start <- pt$start[pt$par_free > 0]

  # Constraints currently only used for measurement invariance
  # in latent variables
  if (constraints@con_logical) {
    # Transform x-vector to "shorter" version
    x_start <- x_start %*% constraints@eq_constraints_Q2
  }





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
    if (any(pt$par[pt$dest == "Sigma_eta" & pt$type == "var"] <= 0)) {
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
  if (fit.model) {
    fit <- nlminb(
      start = x_start,
      objective = objective_function,
      control = list(
        rel.tol = 1e-10 # ,
        # eval.max = 200, # 500
        # iter.max = 150 # 300
      )
    )
  } else {
    fit <- list(
      par = pt$par[pt$par_free > 0],
      convergence = 0
    )
  }




  # Save parameters back to partable
  if (constraints@con_logical & fit.model) {
    pt$par[pt$par_free > 0L] <- as.numeric(
      constraints@eq_constraints_Q2 %*% fit$par
    )
  } else {
    pt$par[pt$par_free > 0L] <- fit$par
  }

  # Compute matrices for checking
  modellist <- creg_modellist(
    pt = pt,
    datalist = datalist,
    gh_grid = gh_grid,
    family = family,
    input = input
  )

  # End of Timer
  if (!silent) {
    time_diff <- Sys.time() - time_start
    units(time_diff) <- "secs"
    cat("done. Took:", round(time_diff, 1), "s\n")
  }

  # Save fit and partable back to object and return object
  object@fit <- new("creg_fit",
    fit = fit,
    objective = objective_function,
    modellist = modellist
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
  gh_grid <- modellist$gh_grid
  cfa <- modellist$cfa

  # if (length(gh_grid$W) > 0) {
  #   # browser()
  # }

  obj_outgroup <- sum(dpois(n_cell, exp(kappas), log = TRUE))

  obj_ingroups <- mapply(function(data, modellist_g) {
    beta <- modellist_g$beta
    Beta <- modellist_g$Beta
    gamma <- modellist_g$gamma
    Gamma <- modellist_g$Gamma
    Omega <- modellist_g$Omega
    overdis <- modellist_g$overdis
    nu <- modellist_g$nu
    Lambda <- modellist_g$Lambda
    Theta <- modellist_g$Theta
    mu_eta <- modellist_g$mu_eta
    Sigma_eta <- modellist_g$Sigma_eta
    mu_z <- modellist_g$mu_z
    Sigma_z <- modellist_g$Sigma_z
    Sigma_z_lv <- modellist_g$Sigma_z_lv
    fixed_z <- modellist_g$fixed_z
    dims <- modellist_g$dims
    N_g <- dims[1]

    if (any(!is.na(Sigma_z))) {
      if (any(diag(solve(Sigma_z)) <= 0)) {
        return(-Inf)
      }
    }
    if (any(!is.na(Theta))) {
      if (any(diag(Theta) <= 0)) {
        return(-Inf)
      }
    }

    if (family == "nbinom") {
      if (overdis <= 0) {
        return(-Inf)
      }
    }

    obj_i <- compute_groupcond_logl(
      y = data$y, w = data$w,
      z = data$z, N = N_g,
      beta = beta, Beta = Beta, gamma = gamma, Gamma = Gamma, Omega = Omega,
      overdis = overdis,
      nu = nu, Lambda = Lambda, Theta = Theta,
      mu_eta = mu_eta, Sigma_eta = Sigma_eta,
      fixeta = gh_grid$X, ghweight = gh_grid$W,
      mu_z = mu_z, Sigma_z = Sigma_z, Sigma_z_lv = Sigma_z_lv,
      fixed_z = fixed_z, cfa = cfa, cores = 1
    )



    if (is.na(obj_i)) {
      return(-Inf)
    } else if (is.infinite(obj_i)) {
      return(-Inf)
    }
    return(obj_i)
  }, data = datalist, modellist_g = modellist$modellist_g, SIMPLIFY = TRUE)

  obj <- -(obj_outgroup + sum(obj_ingroups)) / sum(n_cell)
  if (is.na(obj)) {
    return(+Inf)
  }
  # cat(obj, "\n")
  # if (is.infinite(obj) & sign(obj) == -1) {
  #   browser()
  # }
  return(obj)
}
