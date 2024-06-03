#' Transform Partable to Modellist
#'
#' Modellist contains mu, sigmaw, sigmaz for each group
#'
#' @param pt Parameter table
#' @param datalist Data object (including the data list)
#' @param gh_grid Gauss-Hermite grid
#' @param family Poisson or negative binomial
#' @param input a lavacreg input object
#'
#' @noRd
creg_modellist <- function(pt, datalist, gh_grid, family, input) {
    # Information passed to mapply
    pt_list <- split(pt$par, pt$group)
    n_cell <- input@n_cell
    groupw <- pt$par[pt$dest == "groupw"] # just used for list afterwards

    # Additional information required for actual modellist-function
    modellist_info <- list()
    modellist_info$no_lv <- input@no_lv
    modellist_info$no_w <- input@no_w
    modellist_info$no_z <- input@no_z
    modellist_info$no_int_z <- input@no_int_z
    modellist_info$no_int_lv <- input@no_int_lv
    modellist_info$no_int_z_lv <- input@no_int_z_lv
    modellist_info$no_integration_points <- length(gh_grid$W)

    # row indicators for subsetting
    pt_tmp <- pt[pt$group == 1, ]

    # Regression model
    modellist_info$ind_pt_beta <- grep("beta", pt_tmp$dest)
    modellist_info$ind_pt_Beta <- grep("Beta", pt_tmp$dest)
    modellist_info$ind_pt_gamma <- grep("gamma", pt_tmp$dest)
    modellist_info$ind_pt_Gamma <- grep("Gamma", pt_tmp$dest)
    modellist_info$ind_pt_Omega <- grep("Omega", pt_tmp$dest)
    modellist_info$ind_pt_overdis <- grep("overdis", pt_tmp$dest)

    # Measurement model
    modellist_info$ind_pt_mm_nu <- grep("nu", pt_tmp$dest)
    modellist_info$ind_pt_mm_lambda <- grep("Lambda", pt_tmp$dest)
    modellist_info$ind_pt_sigmaw <- grep("Theta", pt_tmp$dest)

    # LVs
    modellist_info$ind_pt_lv_mu <- grep("mu_eta", pt_tmp$dest)
    modellist_info$ind_pt_lv_var <- intersect(
        grep("Sigma_eta", pt_tmp$dest),
        grep("var", pt_tmp$type)
    )
    modellist_info$ind_pt_lv_cov <- intersect(
        grep("Sigma_eta", pt_tmp$dest),
        grep("cov", pt_tmp$type)
    )

    # Stochastic observed variables
    modellist_info$fixed_z <- FALSE
    if (!is.null(input@creg_options$fixed_z)) {
        modellist_info$fixed_z <- input@creg_options$fixed_z
    }
    modellist_info$ind_pt_z_mu <- grep("mu_z", pt_tmp$dest)
    modellist_info$ind_pt_z_var <- intersect(
        grep("Sigma_z", pt_tmp$dest),
        grep("var", pt_tmp$type)
    )
    modellist_info$ind_pt_z_cov <- intersect(
        grep("Sigma_z", pt_tmp$dest),
        grep("cov", pt_tmp$type)
    )

    modellist_info$ind_pt_z_lv_cov <- grep("Sigma_z_lv", pt_tmp$dest)



    modellist_g <- mapply(modellist_testfun,
        data = datalist,
        pt_g = pt_list,
        N_g = n_cell,
        MoreArgs = list(modellist_info = modellist_info, family = family),
        SIMPLIFY = FALSE
    )

    list(
        groupw = groupw,
        n_cell = n_cell,
        family = family,
        gh_grid = gh_grid,
        modellist_g = modellist_g
    )
}



modellist_testfun <- function(data, pt_g, N_g, modellist_info, family) {
    no_lv <- modellist_info$no_lv
    no_w <- modellist_info$no_w
    no_z <- modellist_info$no_z
    no_int_z <- modellist_info$no_int_z
    no_int_lv <- modellist_info$no_int_lv
    no_int_z_lv <- modellist_info$no_int_z_lv
    no_integration_points <- modellist_info$no_integration_points

    ind_pt_beta <- modellist_info$ind_pt_beta
    ind_pt_Beta <- modellist_info$ind_pt_Beta
    ind_pt_gamma <- modellist_info$ind_pt_gamma
    ind_pt_Gamma <- modellist_info$ind_pt_Gamma
    ind_pt_Omega <- modellist_info$ind_pt_Omega
    ind_pt_overdis <- modellist_info$ind_pt_overdis
    ind_pt_lv_mu <- modellist_info$ind_pt_lv_mu
    ind_pt_lv_var <- modellist_info$ind_pt_lv_var
    ind_pt_lv_cov <- modellist_info$ind_pt_lv_cov
    ind_pt_mm_nu <- modellist_info$ind_pt_mm_nu
    ind_pt_mm_lambda <- modellist_info$ind_pt_mm_lambda
    ind_pt_sigmaw <- modellist_info$ind_pt_sigmaw
    ind_pt_z_mu <- modellist_info$ind_pt_z_mu
    ind_pt_z_var <- modellist_info$ind_pt_z_var
    ind_pt_z_lv_cov <- modellist_info$ind_pt_z_lv_cov
    ind_pt_z_cov <- modellist_info$ind_pt_z_cov

    fixed_z <- modellist_info$fixed_z

    # Extract parameter information
    # 1.1 Regression coefficients
    beta <- pt_g[ind_pt_beta]

    if (no_int_z) {
        Beta_par <- pt_g[ind_pt_Beta]
        Beta <- creg_matrix_vech_reverse(Beta_par, diagonal = FALSE)
        Beta <- Beta * lower.tri(Beta, diag = FALSE)
    } else {
        Beta <- matrix(0, 0, 0)
    }


    if (family == "poisson") {
        overdis <- 0
    } else {
        overdis <- pt_g[ind_pt_overdis]
    }


    if (no_lv) {
        # 1.2 Regression coefficients
        gamma <- pt_g[ind_pt_gamma]

        if (no_int_lv) {
            Gamma_par <- pt_g[ind_pt_Gamma]
            Gamma <- creg_matrix_vech_reverse(Gamma_par, diagonal = FALSE)
            Gamma <- Gamma * lower.tri(Gamma, diag = FALSE)
        } else {
            Gamma <- matrix(0, 0, 0)
        }

        if (no_int_z_lv) {
            Omega_par <- pt_g[ind_pt_Omega]
            Omega <- matrix(Omega_par, ncol = no_lv)
        } else {
            Omega <- matrix(0, 0, 0)
        }

        # 1.3 Latent covariates
        mu_eta <- pt_g[ind_pt_lv_mu]
        lv_var <- pt_g[ind_pt_lv_var]
        lv_cov <- pt_g[ind_pt_lv_cov]

        Sigma_eta <- creg_matrix_vech_reverse(lv_cov, diagonal = FALSE)
        diag(Sigma_eta) <- lv_var

        # 1.4 Measurement model
        nu <- pt_g[ind_pt_mm_nu]
        lambda_vec <- pt_g[ind_pt_mm_lambda]
        Lambda <- matrix(lambda_vec, ncol = no_lv)

        Theta <- matrix(0, no_w, no_w)
        diag(Theta) <- pt_g[ind_pt_sigmaw]
    } else {
        gamma <- numeric()
        Gamma <- matrix(0, 0, 0)
        Omega <- matrix(0, 0, 0)

        mu_eta <- numeric()
        Sigma_eta <- matrix(0, 0, 0)

        nu <- numeric()
        Lambda <- matrix(0, 0, 0)
        Theta <- matrix(0, 0, 0)
    }


    if (no_z & !fixed_z) {
        # 1.4 Manifest covariate
        mu_z <- pt_g[ind_pt_z_mu]
        z_var <- pt_g[ind_pt_z_var]
        z_cov <- pt_g[ind_pt_z_cov]

        Sigma_z <- creg_matrix_vech_reverse(z_cov, diagonal = FALSE)
        diag(Sigma_z) <- z_var
    } else {
        mu_z <- numeric()
        Sigma_z <- matrix(0, 0, 0)
    }

    if (no_z & no_lv & !fixed_z) {
        # 1.5 Manifest-latent covariates covariances
        Sigma_z_lv <- matrix(pt_g[ind_pt_z_lv_cov], ncol = no_lv)
    } else {
        Sigma_z_lv <- matrix(0, 0, 0)
    }


    # 2. dimension information for Cpp
    dims <- c(N_g, no_integration_points, no_z, no_w)

    list(
        beta = beta,
        Beta = Beta,
        gamma = gamma,
        Gamma = Gamma,
        Omega = Omega,
        overdis = overdis,
        nu = nu,
        Lambda = Lambda,
        Theta = Theta,
        mu_eta = mu_eta,
        Sigma_eta = Sigma_eta,
        mu_z = mu_z,
        Sigma_z = Sigma_z,
        Sigma_z_lv = Sigma_z_lv,
        fixed_z = fixed_z,
        dims = dims
    )
}
