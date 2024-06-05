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
    creg_options <- input@creg_options

    fixed_z <- FALSE
    if (!is.null(creg_options$fixed_z)) {
      fixed_z <- creg_options$fixed_z |> as.logical()
    }

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
      res_rows <- pt_g$dest == "beta" |
        pt_g$dest == "gamma" |
        pt_g$dest == "Beta" |
        pt_g$dest == "Gamma" |
        pt_g$dest == "Omega"
      res_cols <- c("rhs", "par", "SE", "zval", "pval")
      res <- pt_g[res_rows, res_cols]
      rownames(res) <- res$rhs
      res <- res[, -1]
      names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
      print(res, digits = 3, print.gap = 3)

      if (family != "poisson") {
        # Print overdispersion parameter if it exists
        cat("\n")
        res_rows <- pt_g$dest == "overdis"
        res_cols <- c("par", "SE", "zval", "pval")
        res <- pt_g[res_rows, res_cols]
        rownames(res) <- "Dispersion"
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits = 3, print.gap = 3)
      }

      if ((no_z & !fixed_z) | no_lv) {
        # Print means and variances of the covariates
        cat("\nMeans:\n")
        res_rows <- pt_g$dest == "mu_z" | pt_g$dest == "mu_eta"
        res_cols <- c("lhs", "par", "SE", "zval", "pval")
        res <- pt_g[res_rows, res_cols]
        rownames(res) <- res$lhs
        res <- res[, -1]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits = 3, print.gap = 3)

        cat("\nVariances:\n")
        res_rows <- pt_g$type == "var" &
          (pt_g$dest == "Sigma_z" | pt_g$dest == "Sigma_eta")
        res_cols <- c("lhs", "par", "SE", "zval", "pval")
        res <- pt_g[res_rows, res_cols]
        rownames(res) <- res$lhs
        res <- res[, -1]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits = 3, print.gap = 3)


        cond_cov_z <- no_z >= 2 & !fixed_z
        cond_cov_lv <- no_lv >= 2
        cond_cov_lv_z <- no_z >= 1 & no_lv >= 1 & !fixed_z


        if (cond_cov_z | cond_cov_lv | cond_cov_lv_z) {
          # Print covariances of covariates
          cat("\nCovariances:\n")
          res_rows <- (pt_g$type == "cov" &
            (pt_g$dest == "Sigma_z" | pt_g$dest == "Sigma_eta")) |
            pt_g$dest == "Sigma_z_lv"
          res_cols <- c("lhs", "op", "rhs", "par", "SE", "zval", "pval")
          res <- pt_g[res_rows, res_cols]
          rownames(res) <- paste(res$lhs, res$op, res$rhs)
          res <- res[, -c(1:3)]
          names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
          print(res, digits = 3, print.gap = 3)
        }
      }


      if (no_lv) {
        # Print measurement model
        cat("\nMeasurement Model:\n")

        cat("Intercepts:\n")
        res_rows <- pt_g$dest == "nu"
        res_cols <- c("lhs", "op", "rhs", "par", "SE", "zval", "pval")
        res <- pt_g[res_rows, res_cols]
        rownames(res) <- paste(res$lhs, res$op, res$rhs)
        res <- res[, -c(1:3)]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits = 3, print.gap = 3)

        cat("\nLoadings:\n")
        res_rows <- pt_g$dest == "Lambda"
        res_cols <- c("lhs", "op", "rhs", "par", "SE", "zval", "pval")
        res <- pt_g[res_rows, res_cols]
        rownames(res) <- paste(res$lhs, res$op, res$rhs)
        res <- res[, -c(1:3)]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits = 3, print.gap = 3)

        cat("\nResidual Variances:\n")
        res_rows <- pt_g$dest == "Theta" & pt_g$type == "var"
        res_cols <- c("lhs", "op", "rhs", "par", "SE", "zval", "pval")
        res <- pt_g[res_rows, res_cols]
        rownames(res) <- paste(res$lhs, res$op, res$rhs)
        res <- res[, -c(1:3)]
        names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
        print(res, digits = 3, print.gap = 3)
      }
    }
  }
)
