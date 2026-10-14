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

      if (family == "nbinom") {
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
        res_rows <- pt_g$dest == "Lambda" &
          (pt_g$par_free > 0L | pt_g$par != 0L)
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
  "summary", signature(object = "lvglm"),
  function(object, digits = max(3L, getOption("digits") - 3L), lavoutput = FALSE) {
    input <- object$input
    N_groups <- input$N_groups
    N_cells <- input$N_cells
    N_coef_g <- input$N_coef_g
    prednames <- input$prednames
    N <- sum(N_cells)
    family <- input$family

    disp_logical <- FALSE
    if (is.na(family$dispersion)) {
      disp_logical <- TRUE
    }

    pt <- data.frame(
      par = object$fit_om$par, SE = NA, zval = NA, pval = NA,
      group = rep(1:N_groups, each = N_coef_g)
    )
    vcov_fit <- object$vcov

    if (input$se & !is.null(vcov_fit)) {
      SE <- sqrt(diag(vcov_fit))

      pt$SE <- SE
      pt$zval <- pt$par / pt$SE
      pt$pval <- 2 * (1 - pnorm(abs(pt$zval)))
    }


    cat("\nCall:\n",
      paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n",
      sep = ""
    )

    for (g in 1:N_groups) {
      pt_g <- pt[pt$group == g, ]

      if (disp_logical) {
        disppar <- pt_g$par[N_coef_g]
        pt_g <- pt_g[-N_coef_g, ]
      } else {
        disppar <- family$dispersion
      }

      if (N_groups > 1L) {
        cat(
          paste0(
            "\n\n--------------------- Group ", g, " ----------------------- \n\n"
          )
        )
      }


      # Print regression coefficients
      cat("Coefficients:\n")
      res_cols <- c("par", "SE", "zval", "pval")
      res <- pt_g[, res_cols]
      rownames(res) <- prednames
      names(res) <- c("Estimate", "SE", "Est./SE", "p-value")
      signif.legend <- ifelse(g == N_groups, TRUE, FALSE)
      printCoefmat(res, digits = digits, signif.legend = signif.legend)

      cat(
        "\n(Dispersion parameter for ", family$family,
        " family taken to be ", format(disppar), ")\n"
      )
    }
    aic <- object$fit_om$objective * N * 2 + input$N_coef * 2
    cat("AIC: ", format(aic, digits = max(4L, digits + 1L)), "\n\n",
      "Number of optimizer iterations: ", object$fit_om$iterations,
      "\n",
      sep = ""
    )

    if (lavoutput) {
      cat(
        "\n\n--------------------- Predictor Model ----------------------- \n\n"
      )
      summary(object$fit_pm)
    }
  }
)




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
# setMethod(
#   "predict", signature(object = "lvglm"),
#   function(object, newdata = NULL, type = c("response", "link")) {
#     type <- match.arg(type)

#     input <- object$input
#     N_groups <- input$N_groups
#     N_coef_g <- input$N_coef_g
#     family <- input$family
#     formula <- input$formula

#     dataext <- object$dataext
#     groupvarext <- object$groupvarext
#     groupvar <- object$groupvar
#     coef <- object$fit_om$par
#     disp_logical <- FALSE


#     if (is.null(newdata)) {
#       X <- object$modelmatrix
#     } else {
#       mf <- model.frame(formula, newdata)
#       mt <- attr(mf, "terms")
#       X <- model.matrix(mt, mf, contrasts = NULL)
#     }


#     if (is.na(family$dispersion)) {
#       disp_logical <- TRUE
#     }
#     if (disp_logical) {
#       coefmat <- matrix(coef, ncol = N_groups)

#       # Delete last row, as it contains the dispersion parameter
#       coefmat <- coefmat[-N_coef_g, ] |> as.matrix(ncol = N_groups)
#     } else {
#       coefmat <- matrix(coef, ncol = N_groups)
#     }


#     linkinv <- family$linkinv
#     etay <- numeric(nrow(dataext))

#     for (g in 1:N_groups) {
#       sel_g_id <- groupvarext == levels(groupvar)[g]
#       etay[sel_g_id] <- (X[sel_g_id, ] %*% coefmat[, g])
#     }

#     if (type == "response") {
#       predy <- linkinv(etay)
#     } else if (type == "link") {
#       predy <- etay
#     }

#     dataext$GHw * predy / sqrt(2)
#     data_t <- data.frame(
#       id = dataext$id,
#       predy = predy
#     )
#     tmp <- setDT(data_t)[, .(predy = sum(predy)), by = id]
#     return(tmp$predy)
#   }
# )
