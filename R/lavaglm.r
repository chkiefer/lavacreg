# Missing features
# - Multigroup models
# - Multiple Outcomes
# - support for multinomial / ordinal regression
# - corrected standard errors

#' @importFrom lavaan lavaanify
#' @importFrom MASS mvrnorm
#' @importFrom lavaan sem
#' @import data.table
lvglm <- function(
  omodel, # Syntax of the outcome model
  pmodel, # Syntax of the predictor model
  data, # Data Frame
  family = gaussian(), # as specified in glm() with a family function
  group = NULL, # if a multi-group model is specified
  verbose = FALSE, # should objective be printed out during optimization
  lavaan_options = list(), # will be passed to lavaan in Step 1
  inttype = "GH",
  Nip = 3,
  se = TRUE
) {
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
    ngroups <- length(levels(groupvar))


    # 1. Lavaanify the outcome model
    pt_om0 <- lavaan::lavaanify(
        omodel,
        meanstructure = TRUE,
        int.ov.free = TRUE,
        fixed.x = TRUE,
        auto.var = TRUE,
        ngroups = ngroups
    )
    # Reduce partable to only include the regression stuff
    # i.e., parameters that are either free or user-fixed
    pt_om <- pt_om0[pt_om0$free > 0 | !is.na(pt_om0$ustart), ]

    # Get names of variables in outcome model
    dvname <- lavaan::lavNames(pt_om, "eqs.y") # which is the outcome
    pnames <- lavaan::lavNames(pt_om, "eqs.x") # which are the predictors
    intnames <- lavaan::lavNames(pt_om, "ov.interaction") # which are the interactions
    # note: the outcome model does not contain information which variables
    # are latent or observed; we derive that from the predictor model

    # Check whether residual variance parameter free
    res_mods <- c("gaussian", "nbinom")
    if (!family$family %in% res_mods) {
        id_res <- pt_om$lhs %in% dvname & pt_om$op == "~~"
        pt_om$free[id_res] <- 0
        pt_om$ustart[id_res] <- 0
    }

    # Adjust numbering of free parameters if necessary
    pt_om$free[pt_om$free > 0] <- 1:length(pt_om$free[pt_om$free > 0])


    # 2. Estimation of the Predictor Model
    default_args <- list(
        model = pmodel,
        data = data,
        group = group,
        meanstructure = TRUE # we rely on extracting means; maybe it would also work without
    )
    fit_pm <- do.call(
        "sem",
        c(default_args, lavaan_options)
    )

    # 3. Get further information from predictor model
    # extract predictor names and types
    # latent variables
    lvnames <- pnames[pnames %in% lavaan::lavNames(fit_pm, "lv")]
    # manifest (stochastic) variables
    ovnames <- pnames[pnames %in% lavaan::lavNames(fit_pm, "ov")]
    # manifest (fixed) variables
    fxnames <- pnames[!pnames %in% c(lvnames, ovnames, intnames)]
    # manifest variables (fixed & stochastic)
    opvnames <- c(ovnames, fxnames)

    # Find out, what kinds of interactions are specified
    inttypes <- sapply(
        strsplit(intnames, ":"),
        function(intvars) {
            if (all(intvars %in% lvnames)) {
                # interaction between LVs
                "Gamma"
            } else if (all(intvars %in% opvnames)) {
                # interaction between OVs
                "Beta"
            } else {
                # Interaction between LVs and OVs
                "Omega"
            }
        }
    )

    # All observed variables in predictor model
    # i.e., predictors as well as indicators
    pm_ovnames <- lavNames(fit_pm, "ov")
    pmnames <- c(lvnames, pm_ovnames)


    # 4. Computation of the person-specific Gauss-Hermite points
    # for numerical integration
    # formula from Wikipedia
    # https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
    # Eaton, Morris L. (1983). Multivariate Statistics: a Vector Space Approach. John Wiley and Sons. pp. 116–117

    # as far as I understand, we use regression factor scores and their variance
    # easiest way should be:
    GH_g <- list(ngroups)
    mubar <- lavPredict(fit_pm, method = "regression", acov = TRUE)
    # Initialize GH grid
    # Pruning not meaningful for one dimension
    # (would remove all but one value)
    no_lv <- length(lvnames)
    if (no_lv > 1) {
        prune <- TRUE
    } else if (no_lv == 1) {
        prune <- FALSE
        mubar <- list(mubar)
    }
    GHgrid <- lavacreg:::creg_init_grid(
        Q = length(lvnames),
        ip = Nip,
        type = inttype,
        prune = prune
    )
    for (g in 1:ngroups) {
        mubar_g <- mubar[[g]]
        if (no_lv == 1L) {
            sigbar_g <- attr(mubar_g, "acov")[[g]]
        } else if (no_lv > 1L) {
            sigbar_g <- attr(mubar, "acov")[[g]]
        }


        # Adapt GH grid to covariance of LVs (same for all units)
        GHgrid_ad <- lavacreg:::creg_adapt_grid(0, sigbar_g, GHgrid)
        colnames(GHgrid_ad$X) <- lvnames

        # Apply individual mean structure (unit-specific)
        GH_g[[g]] <- lapply(1:n_cell[g], function(i) {
            t(mubar_g[i, ] + t(GHgrid_ad$X))
        })
    }

    # 5. Create Group-Specific Model Matrices
    modellist <- list()
    for (g in 1:ngroups) {
        pt_om_g <- pt_om[pt_om$group == g, ]

        pt_om_g$id <- 1:nrow(pt_om_g)
        # beta (intercept and observed predictors)
        id.beta <- pt_om_g$id[
            pt_om_g$lhs == dvname &
                pt_om_g$op %in% c("~1", "~") &
                pt_om_g$rhs %in% c("", opvnames)
        ]
        beta <- pt_om_g$ustart[id.beta]
        attr(beta, "id.free") <- pt_om_g$free[id.beta]

        # gamma (latent variables)
        id.gamma <- pt_om_g$id[
            pt_om_g$lhs == dvname &
                pt_om_g$op %in% c("~") &
                pt_om_g$rhs %in% c("", lvnames)
        ]
        gamma <- pt_om_g$ustart[id.gamma]
        attr(gamma, "id.free") <- pt_om_g$free[id.gamma]

        # Beta (interactions of observed predictors)
        ovintnames <- intnames[inttypes == "Beta"]
        Beta <- Beta_att <- matrix(
            0,
            nrow = length(opvnames), ncol = length(opvnames)
        )
        colnames(Beta) <- rownames(Beta) <-
            colnames(Beta_att) <- rownames(Beta_att) <- opvnames
        id.Beta <- pt_om_g$id[
            pt_om_g$lhs == dvname &
                pt_om_g$op %in% c("~") &
                pt_om_g$rhs %in% c("", ovintnames)
        ]
        for (i in id.Beta) {
            pt_slice <- pt_om_g[i, ]
            tmp0 <- pt_slice$rhs
            tmp1 <- strsplit(tmp0, ":") |> unlist()
            Beta[tmp1[1], tmp1[2]] <- pt_slice$ustart
            Beta_att[tmp1[1], tmp1[2]] <- pt_slice$free
        }
        attr(Beta, "id.free") <- Beta_att

        # Gamma (latent interactions)
        lvintnames <- intnames[inttypes == "Gamma"]
        Gamma <- Gamma_att <- matrix(
            0,
            nrow = length(lvnames), ncol = length(lvnames)
        )
        colnames(Gamma) <- rownames(Gamma) <-
            colnames(Gamma_att) <- rownames(Gamma_att) <- lvnames
        id.Gamma <- pt_om_g$id[
            pt_om_g$lhs == dvname &
                pt_om_g$op %in% c("~") &
                pt_om_g$rhs %in% c("", lvintnames)
        ]
        for (i in id.Gamma) {
            pt_slice <- pt_om_g[i, ]
            tmp0 <- pt_slice$rhs
            tmp1 <- strsplit(tmp0, ":") |> unlist()
            Gamma[tmp1[1], tmp1[2]] <- pt_slice$ustart
            Gamma_att[tmp1[1], tmp1[2]] <- pt_slice$free
        }
        attr(Gamma, "id.free") <- Gamma_att

        # Omega (latent-manifest interactions)
        lvovintnames <- intnames[inttypes == "Omega"]
        Omega <- Omega_att <- matrix(
            0,
            nrow = length(lvnames), ncol = length(opvnames)
        )
        colnames(Omega) <- colnames(Omega_att) <- opvnames
        rownames(Omega) <- rownames(Omega_att) <- lvnames
        id.Omega <- pt_om_g$id[
            pt_om_g$lhs == dvname &
                pt_om_g$op %in% c("~") &
                pt_om_g$rhs %in% c("", lvovintnames)
        ]
        # this should currently only work if LV is mentioned first!
        for (i in id.Omega) {
            pt_slice <- pt_om_g[i, ]
            tmp0 <- pt_slice$rhs
            tmp1 <- strsplit(tmp0, ":") |> unlist()
            Omega[tmp1[1], tmp1[2]] <- pt_slice$ustart
            Omega_att[tmp1[1], tmp1[2]] <- pt_slice$free
        }
        attr(Omega, "id.free") <- Omega_att

        # Residual variance / overdispersion parameter
        id.epsilon <- pt_om_g$id[
            pt_om_g$lhs == dvname &
                pt_om_g$op %in% c("~~") &
                pt_om_g$rhs %in% c(dvname)
        ]
        epsilon <- pt_om_g$ustart[id.epsilon]
        attr(epsilon, "id.free") <- pt_om_g$free[id.epsilon]

        if (ngroups > 1L) {
            group_label <- levels(groupvar)[g]
            group_select <- data[[group]] %in% group_label
            Y_data <- data[group_select, dvname]
            opv_data <- as.matrix(data[group_select, opvnames])
        } else {
            Y_data <- data[, dvname]
            opv_data <- as.matrix(data[, opvnames])
        }


        # Create the whole modellist
        modellist[[g]] <- list(
            beta = beta,
            gamma = gamma,
            Beta = Beta,
            Gamma = Gamma,
            Omega = Omega,
            epsilon = epsilon,
            opv_data = opv_data,
            Y_data = Y_data,
            GH_unit = GH_g[[g]],
            linkinv = family$linkinv,
            cdist = family$family
        )
    }

    # 6. Definition of the objective function
    lglm_objective <- function(x, modellist, verbose = FALSE) {
        obj_g <- numeric(ngroups)
        for (g in 1:ngroups) {
            modellist_g <- modellist[[g]]
            # 6.1 Match current estimates to model matrices
            beta <- modellist_g$beta
            beta[is.na(beta)] <- x[attr(beta, "id.free")]
            gamma <- modellist_g$gamma
            gamma[is.na(gamma)] <- x[attr(gamma, "id.free")]

            Beta <- modellist_g$Beta
            Beta[is.na(Beta)] <- x[attr(Beta, "id.free")]
            Gamma <- modellist_g$Gamma
            Gamma[is.na(Gamma)] <- x[attr(Gamma, "id.free")]
            Omega <- modellist_g$Omega
            Omega[is.na(Omega)] <- x[attr(Omega, "id.free")]

            epsilon <- modellist_g$epsilon
            epsilon[is.na(epsilon)] <- x[attr(epsilon, "id.free")]

            # 6.2 Get some basic information
            # (maybe more useful another place)
            N <- length(modellist_g$Y_data)
            Nip <- nrow(GHgrid$X)

            # 6.3 Extend Y to "replicated" format, for vectorized computations
            Y_data <- modellist_g$Y_data |>
                unlist() |>
                rep(each = Nip)

            GH_unit <- modellist_g$GH_unit


            # 6.4 Compute Linear Predictor and Conditional Mean
            lp_beta <- cbind(modellist_g$opv_data, 1) %*% matrix(beta, ncol = 1) |>
                rep(each = nrow(GH_unit[[1]]))
            #
            lp_gamma <- lapply(
                GH_unit, function(X) X %*% matrix(gamma, ncol = 1)
            ) |>
                unlist()
            #
            lp_Beta <- ((modellist_g$opv_data %*% Beta) * modellist_g$opv_data) |>
                rowSums() |>
                rep(each = nrow(GH_unit[[1]]))
            #
            lp_Gamma <- lapply(
                GH_unit, function(X) rowSums((X %*% Gamma) * X)
            ) |>
                unlist()
            #
            lp_Omega <- lapply(
                1:N, function(i) {
                    GH_unit[[i]] %*% Omega %*% modellist_g$opv_data[i, ]
                }
            ) |> unlist()

            lp <- lp_beta + lp_Beta + lp_gamma + lp_Gamma + lp_Omega
            muy <- modellist_g$linkinv(lp)

            # 6.5 Compute conditional densities
            if (modellist_g$cdist == "poisson") {
                fy <- dpois(Y_data, muy, log = FALSE)
            } else if (modellist_g$cdist == "gaussian") {
                fy <- dnorm(Y_data, muy, sqrt(epsilon))
            } else if (modellist_g$cdist == "nbinom") {
                fy <- dnbinom(Y_data, mu = muy, size = epsilon)
            } else if (modellist_g$cdist == "binomial") {
                fy <- dbinom(Y_data, size = 1, prob = muy)
            }

            # Weigh density by GH weight
            out <- GHgrid$W * fy


            # Sum out within persons using data.table
            data_t <- data.frame(
                id = rep(1:N, each = Nip),
                out = out
            )

            tmp <- setDT(data_t)[, .(x = sum(out)), by = id]

            if (anyNA(tmp)) {
                return(+Inf)
            }

            # 6.6 Logarithmize individual likelihoods and sum
            out <- tmp$x |>
                log() |>
                sum()

            # Scaling and inversion for nlmnib
            obj <- -1 * out
            obj_g[g] <- obj / N
        }

        obj <- sum(obj_g)
        if (verbose) cat("obj = ", obj, "\n")

        return(obj)
    }
    #
    # 7. Set starting values
    npar <- max(pt_om$free) # total number of free parameters
    pt_om$start <- 0 # at the moment, all coefficients are set to zero
    # only residual variances / overdispersion set to 1
    pt_om$start[pt_om$op == "~~"] <- 1


    # 8. Optimization of likelihood
    fit_om <- nlminb(
        start = pt_om$start[pt_om$free > 0],
        objective = lglm_objective,
        modellist = modellist,
        verbose = verbose
    )


    # Write estimates to partable
    pt_om$est[pt_om$free > 0] <- fit_om$par


    if (se) {
        # Compute standard errors
        information <- numDeriv::hessian(
            func = lglm_objective,
            x = fit_om$par,
            modellist = modellist
        )
        vcov_fit <- solve(information) / nrow(data)
        pt_om$se[pt_om$free > 0] <- vcov_fit |>
            diag() |>
            sqrt()
    } else {
        vcov_fit <- matrix()
    }


    # Return some results
    list(pt_om, fit_pm, vcov_fit)
}


nbinom <- function(link = "log") {
    linktemp <- substitute(link)
    if (!is.character(linktemp)) {
        linktemp <- deparse(linktemp)
    }
    okLinks <- c("log", "identity", "sqrt")
    family <- "nbinom"
    if (linktemp %in% okLinks) {
        stats <- make.link(linktemp)
    } else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    } else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name)) {
                linktemp <- stats$name
            }
        } else {
            stop(
                gettextf(
                    "link \"%s\" not available for %s family; available links are %s",
                    linktemp, family, paste(sQuote(okLinks), collapse = ", ")
                ),
                domain = NA
            )
        }
    }
    structure(list(
        family = family, link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv
    ), class = "family")
}
