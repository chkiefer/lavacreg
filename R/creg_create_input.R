#' Create input
#' 
#' Turns the input into an input object
#' 
#' @inheritParams countreg
#' 
#' @importFrom stats as.formula
#' @importFrom stats terms.formula
#' @keywords internal
#' @noRd
creg_create_input <- function(forml, lv, group, data, family, silent, se, cregOptions){
  # Convert formula and extract terms
  forml <- as.formula(forml)
  
  # Extract variable names involved in regression
  vnames <- all.vars(forml)
  dvname <- vnames[1]
  
  if (!is.null(lv)){
    lvnames <- names(lv)
    lvnames <- lvnames[lvnames %in% vnames]
    ovnames <- unname(unlist(lv))
    cvnames <- vnames[!vnames %in% c(lvnames, dvname)]
    # If latent variables are specified, but not added in formula
    if (!length(lvnames)){
      warning("CountReg warning: A latent variable has been specified, but not included in formula.")
      lv <- list()
      lvnames <- character()
      ovnames <- character()
    }
  } else {
    lv <- list()
    lvnames <- character()
    ovnames <- character()
    cvnames <- vnames[!vnames %in% c(lvnames, dvname)]
  }
  
  
  if (is.null(group)){
    group <- character()
  }
  
  if (is.null(cregOptions)){
    cregOptions <- list()
  }
  
  # Do some checking before returning the input object
  # CountReg does not support interactions or higher-order terms
  if (any(attr(terms.formula(forml), "order") != 1)) {
    stop("CountReg Error: Please do not use higher-order terms in your formula.")
  }
  
  # TODO: Check if dv is count
  if (!is.count(data[dvname])){
    stop("CountReg Error: Dependent variable is not a count variable.")
  }
  
  
  res <- new("input",
    forml = forml,
    lvlist = lv,
    vnames = vnames,
    dvname = dvname,
    lvnames = lvnames,
    ovnames = ovnames,
    cvnames = cvnames,
    groupname = group,
    family = family,
    data = data,
    silent = silent,
    se = se,
    cregOptions = cregOptions)
  
  return(res)
}


#' Check for count variable
#' 
#' Checks if the variable is a count variable
#' 
#' @param x vector to be checked
#' @param tol Tolerance
#' @return Function returns logical value indicating whether x can be considered
#' a count variable or not.
#' 
#' @importFrom stats na.omit
#' @export
is.count <- function(x, tol = .Machine$double.eps^0.5){
    x <- na.omit(x)
    x <- unlist(x)
    tmp0 <- abs(x - as.integer(x)) < tol
    tmp1 <- sign(x) == -1
    
    if (sum(tmp1) > 0){
      FALSE
    } else if (sum(tmp0)/length(x) == 1){
      TRUE
    } else  {FALSE}
  }