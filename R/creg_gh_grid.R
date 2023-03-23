creg_gh_grid <- function(input) {
    no_lv <- input@no_lv
    creg_options <- input@creg_options


    # Initialize empty grid of integration points and
    # define number of integration points per dimension
    init_grid <- list()
    if (is.null(creg_options$intPoints) | !is.integer(creg_options$intPoints)) {
        no_integration_points <- 15L
    } else {
        no_integration_points <- creg_options$intPoints
    }

    # If latent variables specified: compute base integration grid
    # (i.e., multidimensional for standard normal)
    if (no_lv) {
        init_grid <- creg_init_grid(
            Q = no_lv,
            ip = no_integration_points,
            type = "GH"
        )
    }

    return(init_grid)
}

#' Create Gauss-Hermite grid
#'
#' Computes a standard Gauss-Hermite grid for Q dimensions with ip integration points
#'
#' @param Q How many dimensions should the grid have?
#' @param ip Number of integration points per dimension
#' @importFrom fastGHQuad gaussHermiteData
#' @importFrom SparseGrid createSparseGrid
#'
#' @noRd
creg_init_grid <- function(Q = 2, ip = 6, type = "GH") {
    if (type == "GH") {
        x <- fastGHQuad::gaussHermiteData(ip)
        w <- x$w / sqrt(pi)
        x <- x$x * sqrt(2)
        X <- as.matrix(expand.grid(lapply(apply(
            replicate(Q, x), 2, list
        ), unlist)))
        g <- as.matrix(expand.grid(lapply(apply(
            replicate(Q, w),
            2, list
        ), unlist)))
        W <- apply(g, 1, function(x) sum(log(x)))
    } else if (type == "Sparse") {
        x <- SparseGrid::createSparseGrid("KPN", Q, 5L)
        X <- x$nodes
        W <- log(x$weights)
        w <- NULL
    }


    return(invisible(list(X = X, W = W, w = w, type = type)))
}


creg_adapt_grid <- function(mu, Sigma, init_grid, prune = FALSE) {
    X <- init_grid$X
    W <- init_grid$W
    w <- init_grid$w
    type <- init_grid$type
    Q <- length(mu)

    trans <- function(X, Sigma) {
        lambda <- with(eigen(Sigma), {
            if (any(values < 0)) {
                warning("Matrix is not positive definite.")
            }
            if (length(values) > 1) {
                vectors %*% diag(sqrt(values))
            } else {
                vectors * sqrt(values)
            }
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
