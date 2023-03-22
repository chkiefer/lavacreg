#' Add linear constraints
#'
#' Adds linear constraints for group-invariant measurement models
#'
#' @param pt Parameter table
#'
#' @noRd
creg_constraints <- function(object) {
  # What I need:
  # constraint matrix A
  # QR <- qr(t(A))
  # ranK <- QR$rank
  # Q <- qr.Q(QR, complete = TRUE)
  # Q2 <- Q[,-seq_len(ranK), drop = FALSE]
  # x.red <- x %*% Q2
  # x <- Q2 %*% x.red

  pt <- object@partable
  no_lv <- object@input@no_lv
  no_groups <- object@input@no_groups

  # some value indicating whether constraints exist
  # right now only for LVs in multigroup
  con_logical <- no_lv > 0L & no_groups >= 2L

  # CODE ONLY FOR LVs in MULTIGROUP
  if (con_logical) {
    # Binding variables locally to the function
    dest <- par_free <- NULL

    no_par <- max(pt$par_free)
    no_groups <- max(pt$group)
    pt_mm <- subset(pt, dest == "mm" & par_free > 0L)
    no_par_mm <- nrow(pt_mm) / no_groups
    no_con <- no_par_mm * (no_groups - 1L)
    A <- matrix(0, ncol = no_par, nrow = no_con)
    par_split <- split(pt_mm, pt_mm$group)

    id.con <- 1L
    for (i in 1:no_par_mm) {
      for (g in 1:(no_groups - 1)) {
        tmp0 <- par_split[[g]]
        tmp1 <- par_split[[g + 1]]
        A[id.con, tmp0$par_free[i]] <- 1
        A[id.con, tmp1$par_free[i]] <- -1
        id.con <- id.con + 1L
      }
    }

    QR <- qr(t(A))
    ranK <- QR$rank
    Q <- qr.Q(QR, complete = TRUE)
    Q2 <- Q[, -seq_len(ranK), drop = FALSE]
  } else {
    Q2 <- matrix()
    A <- matrix()
  }

  # Return new constraints
  res <- new("constraints",
    con_logical = con_logical,
    eq_constraints_Q2 = Q2,
    con_jac = A
  )
  return(res)
}
