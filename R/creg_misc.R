#' Matrix to Vech reverse
#'
#' Turns a Vech to a symmetric matrix. Taken from the lavaan package
#'
#' @param x Vech to be transformed
#' @param diagonal Is the diagonal of the matrix included in x?
#'
#' @noRd
creg_matrix_vech_reverse <- function(x, diagonal = TRUE) {
    if (diagonal) {
        p <- (sqrt(1 + 8 * length(x)) - 1) / 2
    } else {
        p <- (sqrt(1 + 8 * length(x)) + 1) / 2
    }
    S <- numeric(p * p)
    S[creg_matrix_vech_idx(p, diagonal = diagonal)] <- x
    S[creg_matrix_vechru_idx(p, diagonal = diagonal)] <- x
    attr(S, "dim") <- c(p, p)
    S
}

#' Helper for: Matrix to Vech revers
#'
#' Helper function. Taken from the lavaan package
#'
#' @param n dimension of matrix
#' @param diagonal Is the diagonal of the matrix included in x?
#'
#' @noRd
creg_matrix_vech_idx <- function(n = 1L, diagonal = TRUE) {
    n <- as.integer(n)
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n), n, n, byrow = TRUE)
    if (diagonal) {
        which(ROW >= COL)
    } else {
        which(ROW > COL)
    }
}

#' Helper for Matrix to Vech reverse
#'
#' helper function. Taken from the lavaan package
#'
#' @param n number
#' @param diagonal Is the diagonal of the matrix included in x?
#'
#' @noRd
creg_matrix_vechru_idx <- function(n = 1L, diagonal = TRUE) {
    n <- as.integer(n)
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n), n, n, byrow = TRUE)
    tmp <- matrix(seq_len(n * n), n, n, byrow = TRUE)
    if (diagonal) {
        tmp[ROW >= COL]
    } else {
        tmp[ROW > COL]
    }
}
