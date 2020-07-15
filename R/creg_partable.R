#' Create An Initial Partable for Count Regression Model
#' 
#'  The partable serves three purposes:
#'  1. Translating the input into the required parameters
#'  2. The initial \code{par} column contains starting values for each free parameter
#'  3. Holds the values of each iteration within the fitting process and connects them to model
#'  
#'  @param object An object of class \code{countReg}, which contains input and datalist at this point
#'  
#'  @export
creg_partable <- function(object) {
    no_groups <- object@dataobj@no_groups
    no_lv <- object@dataobj@no_lv
    no_w <- object@dataobj@no_w
    no_z <- object@dataobj@no_z
    lv <- object@input@lvlist
    family <- object@input@family
    
    no_cov <- no_lv + no_z
    no_z_lv_covariance <- no_z * no_lv
    no_z_covariance <- no_z * (no_z - 1)/2
    no_lv_covariance <- no_lv * (no_lv - 1)/2
    
    # Parts of the partable
    dest <- NULL
    type <- NULL
    group <- NULL
    par_free <- NULL
    par <- NULL
    
    # IDs and helpers
    par_free_id <- 0L
    
    # create partable
    for (i in 1:no_groups) {
        # Group weights
        dest <- c(dest, "groupw")
        type <- c(type, NA)
        group <- c(group, i)
        par_free_id <- par_free_id + 1L
        par_free <- c(par_free, par_free_id)
        par <- c(par, 0)
        
        # Regression coefficients
        for (j in 0:no_cov) {
            dest <- c(dest, "regcoef")
            type <- c(type, NA)
            group <- c(group, i)
            par_free_id <- par_free_id + 1L
            par_free <- c(par_free, par_free_id)
            par <- c(par, 0.1)
        }
        
        # Measurement model
        no_ind_before <- 0L
        no_ind_after <- no_w
        
        if (no_lv){
            for (j in 1:no_lv) {
                no_ind <- length(lv[[j]])
                no_ind_after <- no_ind_after - no_ind
                
                if (no_ind_before) {
                    for (l in 1:no_ind_before) {
                        dest <- c(dest, "mm")
                        type <- c(type, "lambda")
                        group <- c(group, i)
                        par_free <- c(par_free, 0)
                        par <- c(par, 0)
                    }
                }
                
                for (k in 1:no_ind) {
                    # Intercepts nu
                    dest <- c(dest, "mm")
                    type <- c(type, "nu")
                    group <- c(group, i)
                    
                    # Fix first intercept to zero
                    if (k != 1) {
                        par_free_id <- par_free_id + 1L
                        par_free <- c(par_free, par_free_id)
                        par <- c(par, 0)
                    } else {
                        par_free <- c(par_free, 0)
                        par <- c(par, 0)
                    }
                    
                    # Factor loadings lambda
                    dest <- c(dest, "mm")
                    type <- c(type, "lambda")
                    group <- c(group, i)
                    
                    # Fix first loading to one
                    if (k != 1) {
                        par_free_id <- par_free_id + 1L
                        par_free <- c(par_free, par_free_id)
                        par <- c(par, 1)
                    } else {
                        par_free <- c(par_free, 0)
                        par <- c(par, 1)
                    }
                }
                
                if (no_ind_after) {
                    for (l in 1:no_ind_after) {
                        dest <- c(dest, "mm")
                        type <- c(type, "lambda")
                        group <- c(group, i)
                        par_free <- c(par_free, 0)
                        par <- c(par, 0)
                    }
                }
                
                no_ind_before <- no_ind_before + no_ind
                
            }
        }
        
        # Means and variances of latent variables
        if (no_lv){
            for (j in 1:no_lv) {
                dest <- c(dest, "lv_grid")
                type <- c(type, "mean")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0)
                
                dest <- c(dest, "lv_grid")
                type <- c(type, "var")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 1)
            }
        }
        
        if (no_lv_covariance){
            for (j in 1:no_lv_covariance) {
                dest <- c(dest, "lv_grid")
                type <- c(type, "cov")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0)
            }
        }
        
        
        # Means and variances of manifest covariates
        if (no_z){
            for (j in 1:no_z) {
                dest <- c(dest, "z")
                type <- c(type, "mean")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0)
                
                dest <- c(dest, "z")
                type <- c(type, "var")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 1)
            }
        }
        
        if (no_z_covariance){
            for (j in 1:no_z_covariance) {
                dest <- c(dest, "z")
                type <- c(type, "cov")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0)
            }
        }
        
        
        # Covariances between manifest and latent covariates
        if (no_z_lv_covariance){
            for (j in 1:no_z_lv_covariance) {
                dest <- c(dest, "z")
                type <- c(type, "cov_z_lv")
                group <- c(group, i)
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 0)
            }
        }
        
        
        # Overdispersion and measurement error variance
        if (no_w){
            for (j in 0:no_w) {
                dest <- c(dest, "sigmaw")
                type <- c(type, ifelse(j, "veps", "size"))
                group <- c(group, i)
                if (j) {
                    par_free_id <- par_free_id + 1L
                    par_free <- c(par_free, par_free_id)
                } else if (family == "poisson") {
                    par_free <- c(par_free, 0)
                } else {
                    par_free_id <- par_free_id + 1L
                    par_free <- c(par_free, par_free_id)
                }
                
                par <- c(par, ifelse(j, 1, ifelse(family == "poisson", 0, 1)))
            }
        } else {
            dest <- c(dest, "sigmaw")
            type <- c(type,  "size")
            group <- c(group, i)
            if (family == "poisson"){
                par_free <- c(par_free, 0)
                par <- c(par, 0) 
            } else {
                par_free_id <- par_free_id + 1L
                par_free <- c(par_free, par_free_id)
                par <- c(par, 1) 
            }
                
        }
        
    }
    pt <- data.frame(dest, type, group, par_free, par)
    return(pt)
}
