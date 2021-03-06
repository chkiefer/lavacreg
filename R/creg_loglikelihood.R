#' lavacreg log-likelihood function
#' 
#' Computes the likelihood function using the C++ function for group-conditional likelihoods
#' 
#' @param datalist Datalist
#' @param modellist Modellist
#'  
#' @importFrom stats dpois
#' @keywords internal
#' @noRd
creg_loglikelihood_function <- function(datalist, modellist) {
  # browser()
    kappas <- modellist$groupw
    n_cell <- modellist$n_cell
    no_groups <- length(kappas)
    family <- modellist$family
    
    obj.group <- sum(dpois(n_cell, exp(kappas), log = TRUE))
    
    # obj.ingroups <- creg_group_logl_cpp(datalist, modellist$modellist_g)
    
    obj.ingroups <- mapply(function(data, modellist_g) {
        muy <- modellist_g$muy
        sigmayw <- modellist_g$sigmayw
        muwz <- modellist_g$muwz
        sigmaz <- modellist_g$sigmaz
        ghweight <- modellist_g$ghweight
        detvarz <- modellist_g$detvarz
        dims <- modellist_g$dims

        if (any(!is.na(sigmaz))){
          if (any(diag(solve(sigmaz)) <= 0)) return(-Inf)
        }
        if (any(sigmayw[-1] <= 0)) return(-Inf)
        if (family == "nbinom" & sigmayw[1] <= 0) return(-Inf)

        obj.i <- compute_groupcond_logl(x = data, muy = muy, sigmayw = sigmayw, muwz = muwz,
                                      sigmaz = sigmaz, ghweight = ghweight, detvarz = detvarz,
                                      dims = dims)
        if(is.na(obj.i)) return(-Inf)
        return(obj.i)
    }, data = datalist, modellist_g = modellist$modellist_g, SIMPLIFY = TRUE)
    
    obj <- -(obj.group + sum(obj.ingroups))/sum(n_cell)
    if (is.na(obj)) return(+Inf)
    return(obj)
    
}


#' Fit the lavacreg model
#' 
#' A wrapper for starting values, optimizing loglik and computation of standard errors
#' 
#'  @param object A lavacreg object
#'  
#' @importFrom stats nlminb
#' @importFrom pracma hessian
#' @keywords internal
#' @noRd
creg_fit_model <- function(object) {
    silent <- object@input@silent
    se <- object@input@se
  
  
    pt <- creg_partable(object)
    
    if (object@dataobj@no_lv > 0L){
      if (!silent){
        cat("Computing starting values...")
        time_start <- Sys.time()
      } 
      
      x.start <- creg_starts_lv(object, pt)
      
      if (!silent){
        time_diff <- Sys.time() - time_start
        units(time_diff) <- "secs"
        cat("done. Took:", round(time_diff,1), "s\n")
      } 
    } else {
      x.start <- matrix(pt$par[pt$par_free > 0L], nrow = 1)
    }
    
    # Constraints currently only used for measurement invariance in latent variables
    if (object@dataobj@no_lv > 0L & object@dataobj@no_groups >= 2L){
      tmp <- creg_constraints(pt)
      object@dataobj@eq_constraints_Q2 <- tmp$Q2
      object@dataobj@con_jac <- tmp$con_jac
      x.start <- x.start %*% object@dataobj@eq_constraints_Q2
    }
    
    dataobj <- object@dataobj
    datalist <- object@dataobj@datalist
    family <- object@input@family

    
    objective_function <- function(x) {
      if (anyNA(x)) return(+Inf)
      x <- matrix(x, ncol = 1)
      
      if (object@dataobj@no_lv > 0L & object@dataobj@no_groups >= 2L){
        x <- as.numeric(dataobj@eq_constraints_Q2 %*% x)
      }    
      
      pt$par[pt$par_free > 0L] <- x
      
      if (any(pt$par[pt$dest == "lv_grid" & pt$type == "var"] <= 0)) return(+Inf)
        
      modellist <- creg_modellist(pt = pt, dataobj = dataobj, family = family)
        
      obj <- creg_loglikelihood_function(datalist, modellist)
      # cat(t(round(obj, 3)),"\r")
      return(obj)
    }
    
    if (!silent){
      cat("Fitting the model...")
      time_start <- Sys.time()
    } 
    
    fit <- nlminb(x.start, objective_function, 
                  control = list(rel.tol = 1e-6,
                                 eval.max = 500,
                                 iter.max = 300))
    
    # L-BFGS-S option
    # fit <- optim
    
    if (!silent){
      time_diff <- Sys.time() - time_start
      units(time_diff) <- "secs"
      cat("done. Took:", round(time_diff,1), "s\n")
    } 
    
    if (!fit$convergence & se){
      if (!silent){
        cat("Computing standard errors...")
        time_start <- Sys.time()
      } 
      information <- hessian(objective_function, fit$par)
      eigvals <- eigen(information, symmetric = TRUE,
                       only.values = TRUE)$values
      if(any(eigvals < -1 * .Machine$double.eps^(3/4))) {
        warning("lavacreg WARNING: information matrix is not positive definite; the model may not be identified")
      }
      vcov_fit <- try(solve(information)/sum(object@dataobj@n_cell), silent = TRUE)
      if (!silent){
        time_diff <-  Sys.time() - time_start
        units(time_diff) <- "secs"
        cat("done. Took:", round(time_diff,1), "s\n")
      } 
    } else if (fit$convergence & se){
        vcov_fit <- NULL
        warning("lavacreg warning: Estimation did not converge. Standard errors are not computed.")
    } else (
      vcov_fit <- NULL
    )
    
    if (object@dataobj@no_lv & object@dataobj@no_groups >= 2){
      pt$par[pt$par_free > 0L] <- as.numeric(dataobj@eq_constraints_Q2 %*% fit$par)
      if (!is.null(vcov_fit))  vcov_fit <- dataobj@eq_constraints_Q2 %*% vcov_fit %*% t(dataobj@eq_constraints_Q2)
    } else {
      pt$par[pt$par_free > 0L] <- fit$par
    }
    
    # print(pt)
    object@fit <- list(fit = fit,
                vcov_fit = vcov_fit,
                pt = pt)
    return(object)
}
