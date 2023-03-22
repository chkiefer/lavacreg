#' @title Create a lavacreg datalist object based on the input
#'
#' @description  This function creates the datalist object required for
#' the further computations. It serves the following purposes:
#' 1. Generate model matrix (df with only relevant variables)
#' 2. Split data into groupwise matrices
#'
#' @param input a lavacreg input object
#'
#' @noRd
creg_datalist <- function(input) {
  dvname <- input@dvname
  ovnames <- input@ovnames
  cvnames <- input@cvnames

  groupvar <- input@groupvar
  data <- input@data

  # Select observed variables for model
  model_vars <- c(dvname, ovnames, cvnames)
  model_matrix <- data[model_vars]

  # Split data into group-wise datasets
  datalist <- lapply(split(model_matrix, groupvar), as.matrix)

  return(datalist)
}



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
