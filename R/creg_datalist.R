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
  datalist <- lapply(split(model_matrix, groupvar), function(df_g) {
    y <- df_g[dvname] |> unlist()
    if (length(ovnames)) {
      w <- df_g[ovnames] |> as.matrix()
    } else {
      w <- matrix(0, 0, 0)
    }

    if (length(cvnames)) {
      z <- df_g[cvnames] |> as.matrix()
    } else {
      z <- matrix(0, 0, 0)
    }

    list(y = y, w = w, z = z)
  })

  return(datalist)
}
