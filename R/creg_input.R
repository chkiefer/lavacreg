#' @title Create lavacreg input object from call
#'
#' @description This function turns the input from the lavacreg function
#' into an input object. It serves three purposes:
#' 1. Save input information for reuse in further steps
#' 2. Compute and extract additional information required for further steps
#' 3. Conduct checks (i.e., misspecification, wrong inputs, etc.)
#'
#' @details The input object should contain
#' (a) the call,
#' (b) all newly computed information,
#' (c) maybe a short information of the checks
#'
#' @inheritParams countreg
#'
#' @importFrom stats as.formula
#' @importFrom stats terms.formula
#' @keywords internal
#' @noRd
creg_create_input <- function(forml,
                              lv,
                              group,
                              data,
                              family,
                              silent,
                              se,
                              creg_options) {
  # Convert formula and extract terms
  forml <- as.formula(forml)

  # Extract variable names involved in regression
  vnames <- all.vars(forml)

  # First name (left-hand side) is dependent variable
  dvname <- vnames[1]

  # If latent variables are specified (or something else is passed to lv)
  if (!is.null(lv)) {
    # Extract latent variable names from list
    lvnames <- names(lv)
    lvnames <- lvnames[lvnames %in% vnames]

    # Names of indicator variables
    ovnames <- unname(unlist(lv))

    # Remaining variable names are treated as manifest covariates
    cvnames <- vnames[!vnames %in% c(lvnames, dvname)]

    # If latent variables are specified, but not added in formula
    if (!length(lvnames)) {
      # Give a warning, but nevertheless continue...
      warning(
        "lavacreg warning: A latent variable has been specified, but not
        included in formula."
      )
      # No relevant latent variables (or indicators) for the model
      lv <- list()
      lvnames <- character()
      ovnames <- character()
    }
    # Else: no latent variable specified
  } else {
    # No relevant latent variables (or indicators) for the model
    lv <- list()
    lvnames <- character()
    ovnames <- character()

    # All covariates are manifest covariates
    cvnames <- vnames[!vnames %in% c(lvnames, dvname)]
  }

  # Get number of each covariate type
  no_lv <- length(lvnames)
  no_w <- length(ovnames)
  no_z <- length(cvnames)

  # If no group variable is specified: group name = ""
  # Get number of groups and respective group sizes
  if (is.null(group)) {
    group <- character()
  }

  # If group variable is specified select corresponding variable
  if (length(group)) {
    groupvar <- as.factor(unlist(data[group]))
  } else {
    # Else group variable contains only 1s (single-group model)
    groupvar <- as.factor(rep(1, nrow(data)))
  }

  # Get group sizes and numbers of variables
  n_cell <- as.integer(table(groupvar))
  no_groups <- length(levels(groupvar))


  # at the moment no specific creg_options are implemented
  if (is.null(creg_options)) {
    creg_options <- list()
  }


  #############
  # CHECKS
  #############

  # lavacreg does not support interactions or higher-order terms
  if (any(attr(terms.formula(forml), "order") != 1)) {
    stop(
      "lavacreg Error: Please do not use higher-order terms in your formula."
    )
  }

  # Check if all variables are in data (or empty "")
  if (!all(
    c(dvname, ovnames, cvnames, group) %in% c(names(data), "")
  )) {
    stop("lavacreg Error: At least one of the specified variables
    is not included within your dataset.")
  }

  # Check if dv is count variable
  if (!is_count(data[dvname])) {
    stop("lavacreg Error: Dependent variable is not a count variable.")
  }

  # Return final input object
  res <- new("input",
    forml = forml,
    lvlist = lv,
    vnames = vnames,
    dvname = dvname,
    lvnames = lvnames,
    ovnames = ovnames,
    cvnames = cvnames,
    groupname = group,
    groupvar = groupvar,
    n_cell = n_cell,
    no_groups = no_groups,
    no_lv = no_lv,
    no_w = no_w,
    no_z = no_z,
    family = family,
    data = data,
    silent = silent,
    se = se,
    creg_options = creg_options
  )

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
is_count <- function(x, tol = .Machine$double.eps^0.5) {
  x <- na.omit(x)
  x <- unlist(x)
  tmp0 <- abs(x - as.integer(x)) < tol
  tmp1 <- sign(x) == -1

  if (sum(tmp1) > 0) {
    val <- FALSE
  } else if (sum(tmp0) / length(x) == 1) {
    val <- TRUE
  } else {
    val <- FALSE
  }
  return(val)
}
