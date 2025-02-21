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
  # ## Setting for testing the SAM function
  # forml <- "dv ~ eta1 + eta2 + z11"
  # data <- example01
  # lv <- list(
  #   eta1 = c("z21", "z22"),
  #   eta2 = c("z41", "z42", "z43")
  # )
  # group <- "treat"
  # family <- "poisson"
  # silent <- FALSE
  # se <- TRUE
  # creg_options <- NULL




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
  # SAM VARIANT!!
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

  ## Estimation of Measurement Model
  fit_mm <- creg_cfa(
    lv = lv,
    data = data,
    group = group,
    silent = silent,
    se = se,
    creg_options = creg_options
  )

  mm_dest <- c("nu", "Lambda", "Theta", "mu_eta", "Sigma_eta")
  id_mm <- fit_mm@partable$dest %in% mm_dest
  id_object <- object@partable$dest %in% mm_dest
  object@partable$par[id_object] <- fit_mm@partable$par[id_mm]
  ids_par_free_comlete <- object@partable$par_free
  object@partable$par_free[id_object] <- 0
  no_par_free_rem <- sum(object@partable$par_free != 0)
  object@partable$par_free[object@partable$par_free != 0] <- c(1:no_par_free_rem)

  # because only allowed for measurement model at the moment!
  con_logical_complete <- object@constraints@con_logical
  object@constraints@con_logical <- FALSE


  # model estimation
  object <- creg_model_estimate(object)

  # Standard error estimation
  object@vcov <- creg_vcov(object)

  if (se) {
    ## so far we have the first step and the uncorrected second step
    ## now to correct the SEs of the second step
    object@partable$par_free <- ids_par_free_comlete
    object@constraints@con_logical <- con_logical_complete
    object2 <- creg_model_estimate(object, fit.model = FALSE)

    if (!silent) {
      cat("Computing standard errors...")
      time_start <- Sys.time()
    }

    tmp_mm <- fit_mm@partable$par_free[id_mm & fit_mm@partable$par_free != 0]
    Sigma11 <- fit_mm@vcov[tmp_mm, tmp_mm]
    Sigma22 <- object@vcov
    # information_joint <- creg_vcov(object2, information_only = TRUE)

    testpar <- object2@partable$par[object2@partable$par_free > 0]
    testobj <- object2@fit@objective

    # hessian(testobj, testpar)[sorty, sorty] |> round(3)

    not_mm_ids <- ids_par_free_comlete[!id_object & ids_par_free_comlete != 0]
    mm_ids <- ids_par_free_comlete[id_object & ids_par_free_comlete != 0]
    no_not_mm_ids <- length(not_mm_ids)
    no_mm_ids <- length(mm_ids)
    sorty <- c(mm_ids, not_mm_ids)
    rev_sorty <- order(sorty)

    superobj <- function(x, x_mm) {
      tmp <- c(x_mm, x)[rev_sorty]
      if (object@constraints@con_logical) {
        # Transform x-vector to "shorter" version
        tmp <- tmp %*% object@constraints@eq_constraints_Q2
      }
      testobj(tmp)
    }
    # testobj(testpar)
    superobj(x = testpar[not_mm_ids], x_mm = testpar[mm_ids])
    gradobj <- function(x, x_reg) {
      grad(superobj, x = x_reg, x_mm = x)
    }
    I21 <- jacobian(gradobj, x = testpar[mm_ids], x_reg = testpar[not_mm_ids])
    I12 <- t(I21)
    I22_inv <- Sigma22 * sum(object@input@n_cell)
    # I_joint_sorted <- information_joint[sorty, sorty]
    # I11 <- I_joint_sorted[1:no_mm_ids, 1:no_mm_ids]
    # I22 <- I_joint_sorted[(no_mm_ids + 1):(no_mm_ids + no_not_mm_ids), (no_mm_ids + 1):(no_mm_ids + no_not_mm_ids)]
    # I21 <- I_joint_sorted[(no_mm_ids + 1):(no_mm_ids + no_not_mm_ids), 1:no_mm_ids]
    # I12 <- I_joint_sorted[1:no_mm_ids, (no_mm_ids + 1):(no_mm_ids + no_not_mm_ids)]
    # I22_inv <- solve(I22)

    Sigma_corrected <- Sigma22 + I22_inv %*% I21 %*% Sigma11 %*% I12 %*% I22_inv

    object@vcov <- lavaan:::lav_matrix_bdiag(Sigma11, Sigma_corrected)[rev_sorty, rev_sorty]
    if (!silent) {
      time_diff <- Sys.time() - time_start
      units(time_diff) <- "secs"
      cat("done. Took:", round(time_diff, 1), "s\n")
    }
  }


  # Return all information back to the user
  return(object)
}
