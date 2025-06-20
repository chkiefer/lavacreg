#' Fitting Count Regression Models with Latent Covariates
#'
#' This function is the main function of the package and can be used to estimate
#' latent variable count regression models in one or multiple group(s).
#'
#' @param forml An object of class \code{\link[stats]{formula}} (or one
#' that can be coerced to that class): a symbolic description of the model to
#' be fitted. The details of model specification are given under Details.
#' @param data A data frame containing all variables specified in \code{forml}
#' and/or indicators of the latent variables specified in \code{lv} (if
#' applicable).
#' @param lv A named list, where names of elements represent the names of the
#'  latent variables and each element consists of a character vector containing
#' variable names of indicators for the respective latent variable, e.g.,
#' \code{list(eta1 = c("z1", "z2", "z3"))}.
#' @param group A group variable. If specified, the regression model specified
#' in \code{forml} is estimated as multi-group model (i.e., within each group).
#' @param family A character indicating the family of the generalized linear
#' model to be estimated. At the moment, \code{"poisson"} (for Poisson
#' regression; default) or \code{"nbinom"} (for negative binomial regression)
#' are available.
#' @param silent Logical. Should informations about the estimation process
#' be suppressed? (Defaults to FALSE)
#' @param se Logical. Should standard errors be computed? Defaults to TRUE.
#' (Can take a while for complex models)
#' @param creg_options optional list of additional options for the estimation
#' procedure
#' @return An object of type \code{lavacreg}. Use \code{summary(object)} to
#' print results containing parameter estimates and their standard errors.
#' @examples
#' fit <- countreg(forml = "dv ~ z11", data = example01, family = "poisson")
#' summary(fit)
#' \donttest{
#' fit <- countreg(
#'   forml = "dv ~ eta1 + z11 + z21",
#'   lv = list(eta1 = c("z41", "z42", "z43")),
#'   group = "treat",
#'   data = example01,
#'   family = "poisson"
#' )
#' summary(fit)
#' }
#'
#' @importFrom methods new
#' @export
countreg <- function(forml,
                     data,
                     lv = NULL,
                     group = NULL,
                     family = "poisson",
                     silent = FALSE,
                     se = TRUE,
                     creg_options = NULL) {
  # Initialize new lavacreg object to store and process information
  object <- new("lavacreg")

  # Create, process and save function input
  # TODO: maybe use matchcall or something similar
  object@input <- creg_input(
    forml = forml,
    lv = lv,
    group = group,
    data = data,
    family = family,
    silent = silent,
    se = se,
    creg_options = creg_options
  )

  start.time <- Sys.time()
  # Create partable
  object@partable <- creg_partable(object@input)

  # Constraints
  object@constraints <- creg_constraints(object)

  # Datalist
  # i.e., split data in group-conditional datasets of dv and covariates
  object@datalist <- creg_datalist(object@input)

  # GH grid
  object@gh_grid <- creg_gh_grid(object@input)

  # Start estimation process
  object@partable <- creg_startvals(object)

  # model estimation
  object <- creg_model_estimate(object)

  # Standard error estimation
  object@vcov <- creg_vcov(object)
  end.time <- Sys.time()

  object@time <- end.time - start.time


  # Return all information back to the user
  return(object)
}




#### Structural-After-Measurement (SAM) variant
### only testing at the moment!!!!
countreg_sam <- function(forml,
                         data,
                         lv,
                         group = NULL,
                         family = "poisson",
                         silent = FALSE,
                         se = TRUE,
                         creg_options = NULL) {
  # Initialize new lavacreg object to store and process information
  object <- new("lavacreg")

  start.time <- Sys.time()
  # Create, process and save function input
  # TODO: maybe use matchcall or something similar
  object@input <- creg_input(
    forml = forml,
    lv = lv,
    group = group,
    data = data,
    family = family,
    silent = silent,
    se = se,
    creg_options = creg_options
  )

  # Make SK factorization required for SAM
  if (!is.null(creg_options$factorization)) {
    if (creg_options$factorization == "KM") {
      warning("lavacreg Warning: SAM is only implemented with SK factorization. Switching to that.")
    }
  }
  creg_options$factorization <- "SK"

  # Set default for number of integration points
  if (is.null(creg_options$intPoints)) {
    creg_options$intPoints <- 5L
  }

  # Create partable
  object@partable <- creg_partable(object@input)

  # Constraints
  object@constraints <- creg_constraints(object)

  # Datalist
  # i.e., split data in group-conditional datasets of dv and covariates
  object@datalist <- creg_datalist(object@input)

  # GH grid
  object@gh_grid <- creg_gh_grid(object@input)

  # Start estimation process
  object@partable <- creg_startvals(object)


  # ADD CHECK FOR LV
  # IF NO LV INVOLVED, CFA GIVES AN ERROR RIGHT NOW


  ## Estimation of Predictor Model
  fit_pm <- creg_cfa(
    lv = lv,
    data = data,
    covariates = object@input@cvnames,
    group = group,
    silent = silent,
    se = se,
    creg_options = creg_options
  )

  # Set predictor model estimates as fixed
  # 1. Define which parameters are fixed
  pm_dest <- c(
    "groupw",
    "nu", "Lambda", "Theta",
    "mu_eta", "Sigma_eta",
    "mu_z", "Sigma_z",
    "Sigma_z_lv"
  )

  # 2. Identify these parameters in PM and whole model
  id_pm <- fit_pm@partable$dest %in% pm_dest
  id_object <- object@partable$dest %in% pm_dest

  # 3. Copy the PM estimates to whole model
  object@partable$par[id_object] <- fit_pm@partable$par[id_pm]

  # 4. Remember which parameters where free before fixing
  ids_par_free_comlete <- object@partable$par_free

  # 5. Fix PM parameters in whole model
  object@partable$par_free[id_object] <- 0

  # 6. Generate new count for Outcome Model parameters
  no_par_free_om <- sum(object@partable$par_free != 0)
  object@partable$par_free[object@partable$par_free != 0] <- c(1:no_par_free_om)


  # Deactivate linear constraints for second step
  # These are currently only implemented in measurement model
  # but remember status in whole model
  con_logical_complete <- object@constraints@con_logical
  object@constraints@con_logical <- FALSE

  # ------------------------------------- #
  # - Step 2 - Outcome Model Estimation - #
  # ------------------------------------- #
  # Estimate model on the remaining free parameters
  object <- creg_model_estimate(object)

  # Standard error estimation
  # unadjusted vcov for OM parameters
  # Note: only if object@input@se is true
  object@vcov <- creg_vcov(object)

  # Write free parameters back to partable
  object@partable$par_free <- ids_par_free_comlete

  if (se) {
    ## so far we have the first step and the uncorrected second step vcovs
    ## now to correct the SEs of the second step

    # Reactivate the equality constraints
    object@constraints@con_logical <- con_logical_complete

    # Get objective function for full model
    object2 <- creg_model_estimate(object, fit.model = FALSE)

    if (!silent) {
      cat("Computing standard errors...")
      time_start <- Sys.time()
    }

    # Extract vcov of PM estimates
    tmp_pm <- fit_pm@partable$par_free[id_pm & fit_pm@partable$par_free != 0]
    Sigma11 <- fit_pm@vcov[tmp_pm, tmp_pm]

    # Extract unadjusted vcov of OM estimates
    Sigma22 <- object@vcov

    # Check if adjustment is required
    seadjust <- TRUE
    if (!is.null(creg_options$seadjust)) {
      seadjust <- creg_options$seadjust
    }


    # Identify which parameters are at which position
    om_ids <- ids_par_free_comlete[!id_object & ids_par_free_comlete != 0]
    pm_ids <- ids_par_free_comlete[id_object & ids_par_free_comlete != 0]

    # Compute order of parameters to resort parameters in full model
    sorty <- c(pm_ids, om_ids)
    rev_sorty <- order(sorty)

    if (seadjust == FALSE) {
      # if no adjustment is required, just return unadjusted vcov
      object@vcov <- lavaan:::lav_matrix_bdiag(Sigma11, Sigma22)[rev_sorty, rev_sorty]
    } else if (seadjust == TRUE) {
      # Construct two functions for estimating only part of the Hessian
      # PM and OM parameters are switched
      # First function takes PM est as free
      # second function takes OM est as free
      obj_om_grad <- function(x, x_pm) {
        # This function has the full objective;
        # x are the om par
        # is called via grad function
        # gives grad of om par for fixed pm pars
        tmp <- c(x_pm, x)[rev_sorty]
        if (object@constraints@con_logical) {
          # Transform x-vector to "shorter" version
          tmp <- tmp %*% object@constraints@eq_constraints_Q2
        }

        fullobj <- object2@fit@objective
        fullobj(tmp)
      }

      #
      obj_pm_grad <- function(x, x_om) {
        # this function is called via jacobian
        # x are the pm par
        # produces off-block-diagonal part of Hessian
        grad(obj_om_grad, x = x_om, x_pm = x)
      }



      # Extract full parameter vector
      fullpar <- object2@partable$par[object2@partable$par_free > 0]

      # Compute I21 as part of the hessian (combination of Jacobian and Grad )
      I21 <- jacobian(obj_pm_grad, x = fullpar[pm_ids], x_om = fullpar[om_ids])
      I12 <- t(I21)

      # I22 is computed from step 2 vcov
      I22_inv <- Sigma22 * (-sum(object@input@n_cell))

      # Compute corrected Sigma for step 2
      Sigma_corrected <- Sigma22 +
        I22_inv %*% I21 %*% Sigma11 %*% I12 %*% I22_inv

      object@vcov <- lavaan:::lav_matrix_bdiag(Sigma11, Sigma_corrected)[rev_sorty, rev_sorty]
    }


    if (!silent) {
      time_diff <- Sys.time() - time_start
      units(time_diff) <- "secs"
      cat("done. Took:", round(time_diff, 1), "s\n")
    }
  }
  end.time <- Sys.time()
  object@time <- end.time - start.time


  # Return all information back to the user
  return(object)
}
