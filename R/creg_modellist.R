#' Transform Partable to Modellist
#'
#' Modellist contains mu, sigmaw, sigmaz for each group
#' 
#' @param pt Parameter table
#' @param dataobj Data object (including the data list)
#' @param family Poisson or negative binomial
#' 
#' @noRd
creg_modellist <- function(pt, dataobj, family) {
    # Information passed to mapply
    datalist <- dataobj@datalist
    pt_list <- split(pt$par, pt$group)
    n_cell <- dataobj@n_cell
    groupw <- pt$par[pt$dest == "groupw"]  # just used for list afterwards
    
    # Additional information required for actual modellist-function
    modellist_info <- list()
    modellist_info$no_lv <- dataobj@no_lv
    modellist_info$no_w <- dataobj@no_w
    modellist_info$no_z <- dataobj@no_z
    modellist_info$init_grid <- dataobj@init_grid
    
    
    # row indicators for subsetting
    pt_tmp <- pt[pt$group==1,]
    modellist_info$ind_pt_regcoef <- grep("regcoef", pt_tmp$dest) 
    modellist_info$ind_pt_lv_mu <- intersect(grep("lv_grid", pt_tmp$dest), grep("mean", pt_tmp$type))
    modellist_info$ind_pt_lv_var <- intersect(grep("lv_grid", pt_tmp$dest), grep("var", pt_tmp$type))
    modellist_info$ind_pt_lv_cov <- intersect(grep("lv_grid", pt_tmp$dest), grep("cov", pt_tmp$type))
    modellist_info$ind_pt_mm_nu <- intersect(grep("mm", pt_tmp$dest), grep("nu", pt_tmp$type))
    modellist_info$ind_pt_mm_lambda <- intersect(grep("mm", pt_tmp$dest), grep("lambda", pt_tmp$type))
    modellist_info$ind_pt_sigmaw <- grep("sigmaw", pt_tmp$dest) 
    modellist_info$ind_pt_z_mu <- intersect(grep("z", pt_tmp$dest), grep("mean", pt_tmp$type))
    modellist_info$ind_pt_z_var <- intersect(grep("z", pt_tmp$dest), grep("var", pt_tmp$type))
    modellist_info$ind_pt_z_lv_cov <- intersect(grep("z", pt_tmp$dest), grep("cov_z_lv", pt_tmp$type))
    modellist_info$ind_pt_z_cov <- setdiff(intersect(grep("z", pt_tmp$dest), grep("cov", pt_tmp$type)),  modellist_info$ind_pt_z_lv_cov)
    
     
    
    modellist_g <- mapply(modellist_testfun, 
        data = datalist, 
        pt_g = pt_list, 
        N_g = n_cell, 
        MoreArgs = list(modellist_info = modellist_info),
        SIMPLIFY = FALSE)
    
    list(groupw = groupw, n_cell = n_cell, family = family, modellist_g = modellist_g)
    
    
}

#' #' Transform Partable to Modellist - old version
#' #'
#' #' Modellist contains mu, sigmaw, sigmaz for each group
#' #' 
#' #' @param pt Parameter table
#' #' @param dataobj Data object (including the data list)
#' #' @param family Poisson or negative binomial
#' #' 
#' #' @noRd

# creg_modellist_old <- function(pt, dataobj, family) {
#     datalist <- dataobj@datalist
#     no_groups <- dataobj@no_groups
#     no_lv <- dataobj@no_lv
#     no_w <- dataobj@no_w
#     no_z <- dataobj@no_z
#     n_cell <- dataobj@n_cell
#     groupw <- pt$par[pt$dest == "groupw"]
#     
#     modellist_g <- lapply(1:no_groups, function(g) {
#         
#         pt_g <- subset(pt, group == g)
#         data <- datalist[[g]]
#         N_g <- n_cell[g]
#         n_var <- no_w + no_z + 1
#         
#         
#         # Extract parameter information 1.1 Regression coefficients
#         regcoef <- unlist(subset(pt_g, dest == "regcoef", select = "par"))
#         
#         
#         if (no_lv){
#             # 1.2 Latent covariates
#             lv_mu <- unlist(subset(pt_g, dest == "lv_grid" & type == "mean", select = "par"))
#             lv_var <- unlist(subset(pt_g, dest == "lv_grid" & type == "var", select = "par"))
#             lv_cov <- unlist(subset(pt_g, dest == "lv_grid" & type == "cov", select = "par"))
#             
#             lv_sigma <- creg_matrix_vech_reverse(lv_cov, diagonal = FALSE)
#             diag(lv_sigma) <- lv_var
#             
#             if (no_lv > 1){
#                 lv_gauss_hermite_grid <- MultiGHQuad::init.quad(Q = no_lv, prior = list(mu = lv_mu, Sigma = lv_sigma), ip = 15L, prune = TRUE)
#             } else if (no_lv == 1) {
#                 # Pruning not meaningful for one dimension (would remove all but one value)
#                 lv_gauss_hermite_grid <- MultiGHQuad::init.quad(Q = no_lv, prior = list(mu = lv_mu, Sigma = lv_sigma), ip = 15L, prune = FALSE) 
#                 if (any(eigen(lv_sigma)$values < 0) ) print(lv_sigma)
#             }
#             
#             lv_gauss_hermite_grid$W <- exp(lv_gauss_hermite_grid$W)
#             no_integration_points <- length(lv_gauss_hermite_grid$W)
#             X <- lv_gauss_hermite_grid$X
#             W <- lv_gauss_hermite_grid$W
#             
#             # 1.3 Measurement model
#             nu <- unlist(subset(pt_g, dest == "mm" & type == "nu", select = "par"))
#             lambda_vec <- unlist(subset(pt_g, dest == "mm" & type == "lambda", select = "par"))
#             lambda_matrix <- matrix(lambda_vec, ncol = no_lv)
#             
#             w_sigma <- unlist(subset(pt_g, dest == "sigmaw", select = par))
#         } else {
#             no_integration_points <- 1
#             w_sigma <- unlist(subset(pt_g, dest == "sigmaw", select = par))
#             W <- 1
#         }
#         
#         
#         if (no_z){
#             # 1.4 Manifest covariate
#             z_mu <- unlist(subset(pt_g, dest == "z" & type == "mean", select = "par"))
#             z_var <- unlist(subset(pt_g, dest == "z" & type == "var", select = "par"))
#             z_cov <- unlist(subset(pt_g, dest == "z" & type == "cov", select = "par"))
#             
#             z_sigma <- creg_matrix_vech_reverse(z_cov, diagonal = FALSE)
#             diag(z_sigma) <- z_var
#         }
#         
#         if (no_z & no_lv){
#             # 1.5 Manifest-latent covariates covariances
#             z_lv_cov <- matrix(unlist(subset(pt_g, dest == "z" & type == "cov_z_lv", select = "par")), ncol = no_lv)
#         }
#         
#         
#         # 2 Conditional and unconditional expectations 
#         # 2.1 Conditional expectation of count outcome
#         # TODO: Still kind of ugly with the data type transformations....
#         muy.nrows <- N_g*no_integration_points
#         muy.ncol <- 1 + no_z
#         ID <- rep(1:N_g, each = no_integration_points)
#         model.matrix <- matrix(NA, nrow = muy.nrows, ncol = muy.ncol)
#         model.matrix[ ,1] <- 1
#         
#         if (no_z){
#             model.matrix[ ,2:(no_z + 1)] <- data[ID, (no_w + 2):n_var]
#         }
#         if (no_lv){
#             model.matrix <- as.matrix(cbind(model.matrix, as.data.frame(X)))
#         }
#         
#         
#         #model.matrix <- cbind(1, data[ID, cvnames], as.data.frame(lv_gauss_hermite_grid$X))
#         #model.matrix <- as.matrix(model.matrix)
#         muy <- exp(matrix(model.matrix %*% regcoef, ncol = N_g))
#         
#         # 2.2 Conditional expectation of indicators W
#         if (no_lv){
#             muw <- t(nu + lambda_matrix %*% t(X))
#         }
#         
#         
#         # 2.3 Conditional expectation and variance Z
#         # Note that the inverse is needed in multivariate density
#         if (no_lv & no_z){
#             muz <- t(z_mu + z_lv_cov %*% solve(lv_sigma) %*% (t(lv_gauss_hermite_grid$X) - lv_mu))
#             sigmaz <- solve(z_sigma - z_lv_cov %*% solve(lv_sigma) %*% t(z_lv_cov))
#             detvarz <- 1/det(sigmaz)  # inverse, because matrix has been inversed
#         } else if (no_z & !no_lv){
#             muz <- z_mu
#             sigmaz <- solve(z_sigma)
#             detvarz <- 1/det(sigmaz)  # inverse, because matrix has been inversed
#         } else if (!no_z){
#             muz <- numeric()
#             sigmaz <- matrix()
#             detvarz <- 1
#         }
#         
#         
#         # 3 Prepare output
#         # 3.1 Combine expectations of W and Z
#         if (no_z & no_w){
#             muwz <- cbind(muw, muz)
#         } else if (no_z & !no_w){
#             muwz <- t(as.matrix(muz))
#         } else if (no_w & !no_z){
#             muwz <- muw
#         } else {
#             muwz <- matrix()
#         }
#         
#         
#         # 3.2 dimension information for Cpp
#         dims <- c(N_g, no_integration_points, no_z, no_w)
#         
#         list(muy = muy, sigmayw = w_sigma, muwz = muwz, sigmaz = sigmaz, ghweight = W, detvarz = detvarz, dims = dims)
#     })
#     
#     list(groupw = groupw, n_cell = n_cell, family = family, modellist_g = modellist_g)
# }


#' Matrix to Vech reverse
#' 
#' Turns a Vech to a symmetric matrix. Taken from the lavaan package
#' 
#' @param x Vech to be transformed
#' @param diagonal Is the diagonal of the matrix included in x?
#' 
#' @noRd
creg_matrix_vech_reverse <- function(x, diagonal = TRUE) {
    if (diagonal) {
        p <- (sqrt(1 + 8 * length(x)) - 1)/2
    } else {
        p <- (sqrt(1 + 8 * length(x)) + 1)/2
    }
    S <- numeric(p * p)
    S[creg_matrix_vech_idx(p, diagonal = diagonal)] <- x
    S[creg_matrix_vechru_idx(p, diagonal = diagonal)] <- x
    attr(S, "dim") <- c(p, p)
    S
}

#' Helper for: Matrix to Vech revers
#' 
#' Helper function. Taken from the lavaan package
#' 
#' @param n dimension of matrix
#' @param diagonal Is the diagonal of the matrix included in x?
#' 
#' @noRd
creg_matrix_vech_idx <- function(n = 1L, diagonal = TRUE) {
    n <- as.integer(n)
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n), n, n, byrow = TRUE)
    if (diagonal) 
        which(ROW >= COL) else which(ROW > COL)
}

#' Helper for Matrix to Vech reverse
#' 
#' helper function. Taken from the lavaan package
#' 
#' @param n number
#' @param diagonal Is the diagonal of the matrix included in x?
#' 
#' @noRd
creg_matrix_vechru_idx <- function(n = 1L, diagonal = TRUE) {
    n <- as.integer(n)
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n), n, n, byrow = TRUE)
    tmp <- matrix(seq_len(n * n), n, n, byrow = TRUE)
    if (diagonal) 
        tmp[ROW >= COL] else tmp[ROW > COL]
}

#' Create Gauss-Hermite grid
#' 
#' Computes a standard Gauss-Hermite grid for Q dimensions with ip integration points
#' 
#' @param Q How many dimensions should the grid have?
#' @param ip Number of integration points per dimension
#' 
#' @noRd
creg_init_grid <- function (Q = 2, ip = 6, type = "GH") 
{
    if (type == "GH"){
        x <- fastGHQuad::gaussHermiteData(ip)
        w <- x$w/sqrt(pi)
        x <- x$x * sqrt(2)
        X <- as.matrix(expand.grid(lapply(apply(replicate(Q, x), 2, list), unlist)))
        g <- as.matrix(expand.grid(lapply(apply(replicate(Q, w), 
                                                2, list), unlist)))
        W <- apply(g, 1, function(x) sum(log(x)))   
    } else if (type == "Sparse"){
        x <- SparseGrid::createSparseGrid("KPN", Q, 5L)
        X <- x$nodes
        W <- log(x$weights)
        w <- NULL
    }
    
    
    return(invisible(list(X = X, W = W, w=w, type=type)))
}

creg_adapt_grid <- function (mu, Sigma, init_grid, prune = FALSE) 
{
    X <- init_grid$X
    W <- init_grid$W
    w <- init_grid$w
    type <- init_grid$type
    Q <- length(mu)
    
    trans <- function(X, Sigma) {
        lambda <- with(eigen(Sigma), {
            if (any(values < 0)) 
                warning("Matrix is not positive definite.")
            if (length(values) > 1) 
                vectors %*% diag(sqrt(values))
            else vectors * sqrt(values)
        })
        t(lambda %*% t(X))
    }
    
    X <- trans(X, Sigma)
    X <- t(t(X) + mu)
    
    if (prune & type == "GH") {
        threshold <- log(min(w)^(Q - 1) * max(w))
        relevant <- W >= threshold
        W <- W[relevant]
        X <- X[relevant, , drop = FALSE]
        
    }
    return(invisible(list(X = X, W = W)))
}


modellist_testfun <- function(data, pt_g, N_g, modellist_info) {
    no_lv <- modellist_info$no_lv
    no_w <- modellist_info$no_w
    no_z <- modellist_info$no_z
    init_grid <- modellist_info$init_grid
    
    ind_pt_regcoef <- modellist_info$ind_pt_regcoef
    ind_pt_lv_mu  <- modellist_info$ind_pt_lv_mu 
    ind_pt_lv_var  <- modellist_info$ind_pt_lv_var 
    ind_pt_lv_cov <- modellist_info$ind_pt_lv_cov
    ind_pt_mm_nu  <- modellist_info$ind_pt_mm_nu 
    ind_pt_mm_lambda  <- modellist_info$ind_pt_mm_lambda 
    ind_pt_sigmaw <- modellist_info$ind_pt_sigmaw
    ind_pt_z_mu <- modellist_info$ind_pt_z_mu 
    ind_pt_z_var <-  modellist_info$ind_pt_z_var 
    ind_pt_z_lv_cov  <- modellist_info$ind_pt_z_lv_cov
    ind_pt_z_cov <- modellist_info$ind_pt_z_cov
    
    n_var <- no_w + no_z + 1
    
    # Extract parameter information 1.1 Regression coefficients
    regcoef <- pt_g[ind_pt_regcoef]
    
    if (no_lv){
        # 1.2 Latent covariates
        lv_mu <- pt_g[ind_pt_lv_mu]
        lv_var <- pt_g[ind_pt_lv_var]
        lv_cov <- pt_g[ind_pt_lv_cov]
        
        lv_sigma <- creg_matrix_vech_reverse(lv_cov, diagonal = FALSE)
        diag(lv_sigma) <- lv_var
        
        if (no_lv > 1){
            # lv_gauss_hermite_grid <- MultiGHQuad::init.quad(Q = no_lv, prior = list(mu = lv_mu, Sigma = lv_sigma), ip = 15L, prune = TRUE)
            lv_gauss_hermite_grid <- creg_adapt_grid(mu = lv_mu, Sigma = lv_sigma, init_grid = init_grid, prune = TRUE)
        } else if (no_lv == 1) {
            # Pruning not meaningful for one dimension (would remove all but one value)
            # lv_gauss_hermite_grid <- MultiGHQuad::init.quad(Q = no_lv, prior = list(mu = lv_mu, Sigma = lv_sigma), ip = 15L, prune = FALSE) 
            lv_gauss_hermite_grid <- creg_adapt_grid(mu = lv_mu, Sigma = lv_sigma, init_grid = init_grid, prune = FALSE)
            if (any(eigen(lv_sigma)$values < 0) ) print(lv_sigma)
        }
        
        lv_gauss_hermite_grid$W <- exp(lv_gauss_hermite_grid$W)
        no_integration_points <- length(lv_gauss_hermite_grid$W)
        X <- lv_gauss_hermite_grid$X
        W <- lv_gauss_hermite_grid$W
        
        # 1.3 Measurement model
        nu <- pt_g[ind_pt_mm_nu]
        lambda_vec <- pt_g[ind_pt_mm_lambda]
        lambda_matrix <- matrix(lambda_vec, ncol = no_lv)
        
        w_sigma <- pt_g[ind_pt_sigmaw]
    } else {
        no_integration_points <- 1
        w_sigma <- pt_g[ind_pt_sigmaw]
        W <- 1
    }
    
    
    if (no_z){
        # 1.4 Manifest covariate
        z_mu <- pt_g[ind_pt_z_mu]
        z_var <- pt_g[ind_pt_z_var]
        z_cov <- pt_g[ind_pt_z_cov]
        
        z_sigma <- creg_matrix_vech_reverse(z_cov, diagonal = FALSE)
        diag(z_sigma) <- z_var
    }
    
    if (no_z & no_lv){
        # 1.5 Manifest-latent covariates covariances
        z_lv_cov <- matrix(pt_g[ind_pt_z_lv_cov], ncol = no_lv)
    }
    
    
    # 2 Conditional and unconditional expectations 
    # 2.1 Conditional expectation of count outcome
    # TODO: Still kind of ugly with the data type transformations....
    muy.nrows <- N_g*no_integration_points
    muy.ncol <- 1 + no_z
    ID <- rep(1:N_g, each = no_integration_points)
    model.matrix <- matrix(NA, nrow = muy.nrows, ncol = muy.ncol)
    model.matrix[ ,1] <- 1
    
    if (no_z){
        model.matrix[ ,2:(no_z + 1)] <- data[ID, (no_w + 2):n_var]
    }
    if (no_lv){
        model.matrix <- as.matrix(cbind(model.matrix, as.data.frame(X)))
    }
    
    
    #model.matrix <- cbind(1, data[ID, cvnames], as.data.frame(lv_gauss_hermite_grid$X))
    #model.matrix <- as.matrix(model.matrix)
    muy <- matrix(model.matrix %*% regcoef, ncol = N_g)
    
    # 2.2 Conditional expectation of indicators W
    if (no_lv){
        muw <- t(nu + lambda_matrix %*% t(X))
    }
    
    
    # 2.3 Conditional expectation and variance Z
    # Note that the inverse is needed in multivariate density
    if (no_lv & no_z){
        muz <- t(z_mu + z_lv_cov %*% solve(lv_sigma) %*% (t(lv_gauss_hermite_grid$X) - lv_mu))
        sigmaz <- solve(z_sigma - z_lv_cov %*% solve(lv_sigma) %*% t(z_lv_cov))
        detvarz <- 1/det(sigmaz)  # inverse, because matrix has been inversed
    } else if (no_z & !no_lv){
        muz <- z_mu
        sigmaz <- solve(z_sigma)
        detvarz <- 1/det(sigmaz)  # inverse, because matrix has been inversed
    } else if (!no_z){
        muz <- numeric()
        sigmaz <- matrix()
        detvarz <- 1
    }
    
    
    # 3 Prepare output
    # 3.1 Combine expectations of W and Z
    if (no_z & no_w){
        muwz <- cbind(muw, muz)
    } else if (no_z & !no_w){
        muwz <- t(as.matrix(muz))
    } else if (no_w & !no_z){
        muwz <- muw
    } else {
        muwz <- matrix()
    }
    
    
    # 3.2 dimension information for Cpp
    dims <- c(N_g, no_integration_points, no_z, no_w)
    
    list(muy = muy, sigmayw = w_sigma, muwz = muwz, sigmaz = sigmaz, ghweight = W, detvarz = detvarz, dims = dims)
}
