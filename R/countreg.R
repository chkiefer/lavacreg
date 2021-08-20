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
  object@input <- creg_create_input(
    forml = forml,
    lv = lv,
    group = group,
    data = data,
    family = family,
    silent = silent,
    se = se,
    creg_options = creg_options
  )

  # Create datalist
  # i.e., split data in group-conditional datasets of dv and covariates
  object@dataobj <- creg_create_datalist(object, data)

  # Create partable
  object@partable <- creg_create_partable(object)

  # Start estimation process
  # TODO: seperate model and standard error estimation
  object <- creg_fit_model(object)

  # Return all information back to the user
  return(object)
}