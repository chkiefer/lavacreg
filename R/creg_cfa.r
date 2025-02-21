#' CFA function as a wrapper for countreg
creg_cfa <- function(
    lv,
    data,
    group = NULL,
    silent = FALSE,
    se = FALSE,
    creg_options = NULL) {
    object <- new("lavacreg")

    lvnames <- names(lv)
    forml <- paste(
        ". ~ ",
        paste(lvnames, collapse = " + ")
    )

    object <- countreg(
        forml = forml,
        data = data,
        lv = lv,
        group = group,
        family = "poisson",
        silent = silent,
        se = se,
        creg_options = creg_options
    )

    return(object)
}
