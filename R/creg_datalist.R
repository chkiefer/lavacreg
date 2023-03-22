#' @title Create a lavacreg datalist object based on the input
#'
#' @description  This function creates the datalist object required for
#' the further computations. It serves the following purposes:
#' 1. Generate model matrix (df with only relevant variables)
#' 2. Split data into groupwise matrices
#' 3. Generate GH grid (if latent variables are specified)
#' # TODO: maybe 1 and 2 suffice, and this would be more of a step
#' # for data wrangling
#'
#' @param object a lavacreg object
#' @param data the dataframe
#'
#' @noRd
creg_create_datalist <- function(object) {
  input <- object@input
  dvname <- input@dvname
  lvnames <- input@lvnames
  ovnames <- input@ovnames
  cvnames <- input@cvnames
  groupvar <- input@groupvar
  no_lv <- input@no_lv
  family <- input@family
  data <- input@data
  creg_options <- input@creg_options


  # Select observed variables for model
  model_vars <- c(dvname, ovnames, cvnames)
  model_matrix <- data[model_vars]

  # Split data into group-wise datasets
  datalist <- lapply(split(model_matrix, groupvar), as.matrix)

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

  # Return new dataobj
  # TODO: placeholder for constraints, maybe not a suitable location
  # TODO: maybe own section for constraints, for increasing complexity
  res <- new("dataobj",
    datalist = datalist,
    eq_constraints_Q2 = matrix(),
    con_jac = matrix(),
    init_grid = init_grid
  )

  return(res)
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
