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
