#' Summary of a lavacreg object
#' 
#' Exports the parameter table with parameter estimates and standard errors 
#' for an estimated latent variable count regression model.
#' 
#' @export
#' @param object lavacreg object
#' @return Function prints the parameter table of an estimated model, which
#' includes the parameter estimates and standard errors.
setMethod("summary", signature(object="lavacreg"),
  function(object) {
    pt <- object@fit$pt
    pt$SE <- NULL
    SE <- sqrt(diag(object@fit$vcov_fit))
    
    pt$SE[as.logical(pt$par_free)] <- SE
    print(pt)
  }
)