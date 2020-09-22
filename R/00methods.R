#' @export
setMethod("summary", signature(object = "lavacreg"),
  function(object, ...) {
    pt <- object@fit$pt
    pt$SE <- NULL
    SE <- sqrt(diag(object@fit$vcov_fit))
    
    pt$SE[as.logical(pt$par_free)] <- SE
    print(pt)
  }
)