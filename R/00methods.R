#' Summary of a lavacreg object
#' 
#' Exports the parameter table with parameter estimates and standard errors 
#' for an estimated latent variable count regression model.
#' 
#' @export
#' @param object lavacreg object
#' @return Function prints the parameter table of an estimated model, which
#' includes the parameter estimates and standard errors.
#' @importFrom stats pnorm
setMethod("summary", signature(object="lavacreg"),
  function(object) {
    input <- object@input
    dataobj <- object@dataobj
    pt <- object@fit$pt
    pt$pval <- pt$zval <- pt$SE <- NA
    vcov_fit <- object@fit$vcov_fit
    
    if (object@input@se & !is.null(vcov_fit)){
      SE <- sqrt(diag(object@fit$vcov_fit))
      
      pt$SE[as.logical(pt$par_free)] <- SE
      pt$zval <- pt$par/pt$SE
      pt$pval <- 2*(1-pnorm(abs(pt$zval)))
    }
    
    for (g in 1:dataobj@no_groups){
      pt_g <- pt[pt$group == g,]
      cat(paste0("\n\n--------------------- Group ",g," --------------------- \n\n"))
      
      # Print regression coefficients
      cat("Regression:\n")
      res <- subset(pt_g, pt_g$dest == "regcoef", select = c("rhs", "par", "SE", "zval", "pval"))
      rownames(res) <- res$rhs
      res <- res[,-1]
      names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
      print(res, digits=3, print.gap=3)
      
      if (input@family != "poisson"){
        # Print overdispersion parameter if it exists
        res <- subset(pt_g, pt_g$type == "size", select = c("par", "SE", "zval", "pval"))
        rownames(res) <- "Dispersion"
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits=3, print.gap=3)
      }
      
      if (dataobj@no_z | dataobj@no_lv){
        # Print means and variances of the covariates
        cat("\nMeans:\n")
        res <- subset(pt_g, pt_g$type == "mean", select = c("lhs", "par", "SE", "zval", "pval"))
        rownames(res) <- res$lhs
        res <- res[,-1]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits=3, print.gap=3)
        
        cat("\nVariances:\n")
        res <- subset(pt_g, pt_g$type == "var", select = c("lhs", "par", "SE", "zval", "pval"))
        rownames(res) <- res$lhs
        res <- res[,-1]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits=3, print.gap=3)
      }
      
      if (dataobj@no_z + dataobj@no_lv >= 2){
        # Print covariances of covariates
        cat("\nCovariances:\n")
        res <- subset(pt_g, pt_g$type == "cov" | pt_g$type == "cov_z_lv", 
          select = c("lhs", "op", "rhs", "par", "SE", "zval", "pval"))
        rownames(res) <- paste(res$lhs, res$op, res$rhs)
        res <- res[,-c(1:3)]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits=3, print.gap=3)
      }
      
      if (dataobj@no_lv){
        # Print measurement model
        cat("\nMeasurement Model:\n")
        res <- subset(pt_g, (pt_g$dest == "mm" | pt_g$type == "veps") & pt_g$par_free > 0, 
          select = c("lhs", "op", "rhs", "par", "SE", "zval", "pval"))
        rownames(res) <- paste(res$lhs, res$op, res$rhs)
        res <- res[,-c(1:3)]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits=3, print.gap=3)
      }
      
    }
    
  }
)