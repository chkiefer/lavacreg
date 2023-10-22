#' Summary of a lavacreg object
#'
#' Exports the parameter table with parameter estimates and standard errors
#' for an estimated latent variable count regression model.
#'
#' @export
#' @param object lavacreg object
#' @return Function prints the parameter table of an estimated model, which
#' includes the parameter estimates and standard errors.
#' @importFrom stats pnorm
setMethod(
  "summary", signature(object = "lavacreg"),
  function(object) {
    input <- object@input
    no_groups <- input@no_groups
    no_z <- input@no_z
    no_lv <- input@no_lv

    family <- input@family
    pt <- object@partable
    pt$pval <- pt$zval <- pt$SE <- NA
    vcov_fit <- object@vcov

    if (input@se & !is.null(vcov_fit)) {
      SE <- sqrt(diag(vcov_fit))

      pt$SE[as.logical(pt$par_free)] <- SE
      pt$zval <- pt$par / pt$SE
      pt$pval <- 2 * (1 - pnorm(abs(pt$zval)))
    }

    for (g in 1:no_groups) {
      pt_g <- pt[pt$group == g, ]
      cat(
        paste0(
          "\n\n--------------------- Group ", g, " ----------------------- \n\n"
        )
      )

      # Print regression coefficients
      cat("Regression:\n")
      res <- subset(
        x = pt_g,
        subset = dest == "beta" | dest == "gamma",
        select = c("rhs", "par", "SE", "zval", "pval")
      )
      rownames(res) <- res$rhs
      res <- res[, -1]
      names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
      print(res, digits = 3, print.gap = 3)

      if (family != "poisson") {
        # Print overdispersion parameter if it exists
        cat("\n")
        res <- subset(
          x = pt_g,
          subset = dest == "overdis",
          select = c("par", "SE", "zval", "pval")
        )
        rownames(res) <- "Dispersion"
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits = 3, print.gap = 3)
      }

      if (no_z | no_lv) {
        # Print means and variances of the covariates
        cat("\nMeans:\n")
        res <- subset(
          x = pt_g,
          subset = dest == "mu_z" | dest == "mu_eta",
          select = c("lhs", "par", "SE", "zval", "pval")
        )
        rownames(res) <- res$lhs
        res <- res[, -1]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits = 3, print.gap = 3)

        cat("\nVariances:\n")
        res <- subset(
          x = pt_g,
          subset = type == "var" & (dest == "Sigma_z" | dest == "Sigma_eta"),
          select = c("lhs", "par", "SE", "zval", "pval")
        )
        rownames(res) <- res$lhs
        res <- res[, -1]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits = 3, print.gap = 3)
      }

      if (no_z + no_lv >= 2) {
        # Print covariances of covariates
        cat("\nCovariances:\n")
        res <- subset(
          x = pt_g,
          subset = (type == "cov" &
            (dest == "Sigma_z" | dest == "Sigma_eta")) | dest == "Sigma_z_lv",
          select = c("lhs", "op", "rhs", "par", "SE", "zval", "pval")
        )
        rownames(res) <- paste(res$lhs, res$op, res$rhs)
        res <- res[, -c(1:3)]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits = 3, print.gap = 3)
      }

      if (no_lv) {
        # Print measurement model
        cat("\nMeasurement Model:\n")

        cat("Intercepts:\n")
        res <- subset(
          x = pt_g,
          subset = dest == "nu",
          select = c("lhs", "op", "rhs", "par", "SE", "zval", "pval")
        )
        rownames(res) <- paste(res$lhs, res$op, res$rhs)
        res <- res[, -c(1:3)]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits = 3, print.gap = 3)

        cat("\nLoadings:\n")
        res <- subset(
          x = pt_g,
          subset = dest == "Lambda",
          select = c("lhs", "op", "rhs", "par", "SE", "zval", "pval")
        )
        rownames(res) <- paste(res$lhs, res$op, res$rhs)
        res <- res[, -c(1:3)]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits = 3, print.gap = 3)

        cat("\nResidual Variances:\n")
        res <- subset(
          x = pt_g,
          subset = dest == "Theta" & type == "var",
          select = c("lhs", "op", "rhs", "par", "SE", "zval", "pval")
        )
        rownames(res) <- paste(res$lhs, res$op, res$rhs)
        res <- res[, -c(1:3)]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits = 3, print.gap = 3)
      }
    }
  }
)
