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
creg_input <- function(forml,
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

  cfa <- FALSE
  if (dvname == ".") {
    cfa <- TRUE
    dvname <- character()
    vnames <- vnames[-1]
  }

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

    if (length(lvnames) != length(names(lv))) {
      # Give a warning, but nevertheless continue...
      warning(
        "lavacreg warning: An additional latent variable has been specified,
         but not included in formula."
      )
      # No relevant latent variables (or indicators) for the model
      lv <- lv[lvnames]
      ovnames <- unname(unlist(lv))
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

  # Interactions
  intnames <- list()
  no_int_z <- 0L
  no_int_lv <- 0L
  no_int_z_lv <- 0L
  forml_terms <- terms.formula(forml)
  int_ind <- attr(forml_terms, "order") == 2

  if (any(int_ind)) {
    int_var <- attr(forml_terms, "factors")[, int_ind] |> as.matrix()
    int_var_logic <- array(as.logical(int_var), dim(int_var))
    int_table <- apply(int_var_logic, 2, function(x) {
      vnames_int <- rownames(int_var)[x]
      if (all(vnames_int %in% cvnames)) {
        return(c("z", vnames_int))
      } else if (all(vnames_int %in% lvnames)) {
        return(c("lv", vnames_int))
      } else {
        return(c("z_lv", vnames_int))
      }
    }) |> t()


    intnames$z <- int_table[int_table[, 1] == "z", 2:3] |> matrix(ncol = 2)
    intnames$lv <- int_table[int_table[, 1] == "lv", 2:3] |> matrix(ncol = 2)
    intnames$z_lv <- int_table[int_table[, 1] == "z_lv", 2:3] |>
      matrix(ncol = 2)

    # Get number of each interaction type
    no_int_z <- nrow(intnames$z)
    no_int_lv <- nrow(intnames$lv)
    no_int_z_lv <- nrow(intnames$z_lv)
  }




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

  # lavacreg does not support higher-order terms
  if (any(attr(terms.formula(forml), "order") > 2)) {
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
  if (cfa) {
    warning("No dependent variable specified; computing a CFA.")
  } else if (!is_count(data[dvname])) {
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
    intnames = intnames,
    groupname = group,
    groupvar = groupvar,
    n_cell = n_cell,
    no_groups = no_groups,
    no_lv = no_lv,
    no_w = no_w,
    no_z = no_z,
    no_int_z = no_int_z,
    no_int_lv = no_int_lv,
    no_int_z_lv = no_int_z_lv,
    family = family,
    data = data,
    silent = silent,
    se = se,
    cfa = cfa,
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
