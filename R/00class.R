#' Input object
#'
#' Takes the lavacreg input
#'
#' @noRd
setClass(
  "input",
  representation(
    forml = "formula",
    lvlist = "list",
    vnames = "character",
    dvname = "character",
    lvnames = "character",
    ovnames = "character",
    cvnames = "character",
    intnames = "list",
    groupname = "character",
    groupvar = "factor",
    n_cell = "integer",
    no_groups = "integer",
    no_lv = "integer",
    no_w = "integer",
    no_z = "integer",
    no_int_z = "integer",
    no_int_lv = "integer",
    no_int_z_lv = "integer",
    family = "character",
    data = "data.frame",
    silent = "logical",
    se = "logical",
    cfa = "logical",
    creg_options = "list"
  )
)

#' Data object
#'
#' Takes the lavacreg data structure
#'
#' @noRd
setClass(
  "dataobj",
  representation(
    datalist = "list",
    eq_constraints_Q2 = "matrix",
    con_jac = "matrix",
    init_grid = "list"
  )
)

#' Constraints object
#'
#' Takes the lavacreg constraints
#'
#' @noRd
setClass(
  "constraints",
  representation(
    con_logical = "logical",
    eq_constraints_Q2 = "matrix",
    con_jac = "matrix"
  )
)

#' Fit object
#'
#' Takes the lavacreg fit information
#'
#' @noRd
setClass(
  "creg_fit",
  representation(
    fit = "list",
    objective = "function",
    modellist = "list"
  )
)


setClassUnion("MatrixOrNULL", c("matrix", "NULL"))

#' lavacreg object
#'
#' The overall object holding all information
#'
#' @noRd
setClass(
  "lavacreg",
  representation(
    input = "input",
    partable = "data.frame",
    constraints = "constraints",
    datalist = "list",
    gh_grid = "list",
    x_start = "matrix",
    fit = "creg_fit",
    vcov = "MatrixOrNULL",
    time = "difftime"
  )
)
