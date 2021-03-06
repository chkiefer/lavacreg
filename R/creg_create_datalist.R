#' Create a datalist object
#' 
#' Creates the dataobject required for computations
#' 
#' @param object a lavacreg object
#' @param data the dataframe
#' 
#' @noRd
creg_create_datalist <- function(object, data){
  input <- object@input
  dvname <- input@dvname
  lvnames <- input@lvnames
  ovnames <- input@ovnames
  cvnames <- input@cvnames
  groupname <- input@groupname
  family <- input@family
  cregOptions <- input@cregOptions
  
  if (length(groupname)){
    groupvar <- as.factor(unlist(data[groupname]))
  } else {
    groupvar <- as.factor(rep(1, nrow(data)))
  }
  
  model.vars <- c(dvname, ovnames, cvnames)
  model.matrix <- data[model.vars]
  datalist <- lapply(split(model.matrix, groupvar), as.matrix)
  
  n_cell <- as.integer(table(groupvar))
  no_groups <- length(levels(groupvar))
  no_lv <- length(lvnames)
  no_w <- length(ovnames)
  no_z <- length(cvnames)
  
  init_grid <- list()
  if (is.null(cregOptions$intPoints) | !is.integer(cregOptions$intPoints)){
    no_integration_points <- 15L
  } else {
    no_integration_points <- cregOptions$intPoints
  }
  
  if (no_lv) init_grid <- creg_init_grid(Q = no_lv, ip = no_integration_points, type = "GH")
  
  
  res <- new("dataobj",
    datalist           = datalist,
    groupvar           = groupvar,
    n_cell = n_cell,
    no_groups = no_groups,
    no_lv = no_lv,
    no_w = no_w,
    no_z = no_z,
    eq_constraints_Q2 = matrix(),
    con_jac = matrix(),
    init_grid = init_grid)
  
  return(res)
}