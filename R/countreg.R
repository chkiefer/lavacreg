#' Fitting Count Regression Models with Latent Covariates
#' 
#' Description
#' 
#' Details
#' 
#' @param forml an object of class for \code{\link[stats]{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under Details.
#' @param lv Definition of the latent variables.
#' @param group A group variable
#' @param data a data frame
#' @param family Poisson or negative binomial
#' @param silent Should informations about the estimation process be suppressed?
#' @param se Should standard errors be computed? (Can take a while for complex models)
#' @return An object of type \code{countreg}.
#' @examples
#' \dontrun{
#' fit <- countreg(forml = 'dv ~ eta1 + z41 + z42', 
#'               lv = list(eta1=c('z11', 'z12')),
#'               group = 'treat',
#'                data = example01,
#'              family = 'poisson')
#' summary(fit)
#' }
#' 
#' 
#' @export
countreg <- function(forml, lv = NULL, group = NULL, data, family, silent = FALSE, se = TRUE) {
  
  object <- new("lavacreg")
  object@input <- creg_create_input(forml, lv, group, data, family, silent, se)
    
  # Create datalist
  # i.e., split data in group-conditional datasets of dv and covariates
  object@dataobj <- creg_create_datalist(object, data)
    
  object <- creg_fit_model(object)
  
  return(object)
}
