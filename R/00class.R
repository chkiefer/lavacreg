#' @export
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
    se = "logical"
    )
)


#' @export
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
    con_jac = "matrix"
  )
)




#' @export
setClass("lavacreg",
  representation(
    input              = "input",
    dataobj           = "dataobj",
    fit = "list"
  )
)