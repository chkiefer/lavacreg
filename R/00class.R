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
    family = "character"
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
    no_z = "integer"
  )
)




#' @export
setClass("countReg",
  representation(
    input              = "input",
    dataobj           = "dataobj",
    fit = "list"
  )
)