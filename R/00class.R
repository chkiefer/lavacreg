#' Input object 
#' 
#' Takes the lavacreg input
#'  
#' @noRd
setClass("input",
  representation(
    forml = "formula",
    lvlist = "list",
    vnames = "character",
    dvname = "character",
    lvnames = "character",
    ovnames = "character",
    cvnames = "character",
    groupname = "character",
    family = "character",
    data = "data.frame",
    silent = "logical",
    se = "logical",
    cregOptions = "list"
    )
)

#' Data object 
#' 
#' Takes the lavacreg data structure
#'  
#' @noRd
setClass("dataobj",
  representation(
    datalist           = "list",
    groupvar           = "factor",
    n_cell = "integer",
    no_groups = "integer",
    no_lv = "integer",
    no_w = "integer",
    no_z = "integer",
    eq_constraints_Q2 = "matrix",
    con_jac = "matrix",
    init_grid = "list"
  )
)

#' lavacreg object 
#' 
#' The overall object holding all information
#'  
#' @noRd
setClass("lavacreg",
  representation(
    input              = "input",
    dataobj           = "dataobj",
    fit = "list"
  )
)