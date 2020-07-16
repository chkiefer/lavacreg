creg_loglikelihood_function <- function(datalist, modellist) {
    kappas <- modellist$groupw
    n_cell <- modellist$n_cell
    no_groups <- length(kappas)
    family <- modellist$family
    
    obj.group <- sum(dpois(n_cell, exp(kappas), log = TRUE))
    
    obj.ingroups <- sapply(1:no_groups, function(g) {
        data <- datalist[[g]]
        modellist_g <- modellist$modellist_g[[g]]
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
        
        obj.i <- CountReg::compute_groupcond_logl(x = data, muy = muy, sigmayw = sigmayw, muwz = muwz, 
                                      sigmaz = sigmaz, ghweight = ghweight, detvarz = detvarz, dims = dims)
        return(obj.i)
    })
    
    obj <- -(obj.group + sum(obj.ingroups))/sum(n_cell)
    return(obj)
    
}



creg_fit_model <- function(object) {
    
    pt <- creg_partable(object)
    dataobj <- object@dataobj
    datalist <- object@dataobj@datalist
    family <- object@input@family
    x.start <- as.numeric(pt$par[pt$par_free != 0])
    
    objective_function <- function(x) {
        if (anyNA(x)) return(+Inf)
        pt$par[as.logical(pt$par_free)] <- x
        
        if (any(pt$par[pt$dest == "lv_grid" & pt$type == "var"] <= 0)) return(+Inf)
        
        modellist <- creg_modellist(pt = pt, dataobj = dataobj, family = family)
        
        substart <- Sys.time()
        obj <- creg_loglikelihood_function(datalist, modellist)
        cpp_time <- round(Sys.time() - substart,3)
        
        # duration <- (Sys.time() - start.time)
        # units(duration) <- "secs"
        # cat("objective:", obj, "time:", round(duration, 1), "cpp time:", cpp_time,"\n")
        return(obj)
    }
    
    cat("Fitting the model.\n")
    start.time <- Sys.time()
    fit <- nlminb(x.start, objective_function, 
                  control = list(rel.tol = 1e-6,
                                 eval.max = 500,
                                 iter.max = 300))
    if (!fit$convergence){
        cat("Computing standard errors.\n")
      vcov_fit <- solve(pracma::hessian(objective_function, fit$par))/sum(object@dataobj@n_cell)
    } else {
        vcov_fit <- NULL
        warning("CountReg warning: Estimation did not converge. Standard errors are not computed.")
    }
    
    pt$par[as.logical(pt$par_free)] <- fit$par
    # print(pt)
    res <- list(fit = fit,
                vcov_fit = vcov_fit,
                pt = pt)
    return(res)
}
