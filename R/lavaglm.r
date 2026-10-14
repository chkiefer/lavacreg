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
  formula, # Formula for the outcome model
  pmodel = NULL, # Syntax of the predictor model
  data, # Data Frame
  family = gaussian(), # as specified in glm() with a family function
  group = NULL, # if a multi-group model is specified
  lavaan_options = list(), # will be passed to lavaan in Step 1
  inttype = "GH",
  Nip = 3,
  se = TRUE,
  verbose = FALSE # should objective be printed out during optimization
) {
    #################################
    ## adapted from glm function
    ############################
    cl <- match.call()

    N <- nrow(data)

    tt <- terms(formula)
    vnames <- attr(tt, "variables") |> as.character()
    vnames <- vnames[-1]


    creg_dpois <- function(Y, mu, disppar) {
        dpois(x = Y, lambda = mu, log = FALSE)
    }

    creg_dnorm <- function(Y, mu, disppar) {
        dnorm(x = Y, mean = mu, sd = sqrt(disppar), log = FALSE)
    }

    creg_dbinom <- function(Y, mu, disppar) {
        dbinom(x = Y, size = 1, prob = mu, log = FALSE)
    }

    creg_dnbinom <- function(Y, mu, disppar) {
        dnbinom(x = Y, mu = mu, size = disppar, log = FALSE)
    }


    family$dfun <- switch(family$family,
        poisson = creg_dpois,
        gaussian = creg_dnorm,
        binomial = creg_dbinom,
        nbinom = creg_dnbinom
    )

    # 6.5 Compute conditional densities
    # if (modellist_g$cdist == "poisson") {
    #     fy <- dpois(Y_data, muy, log = FALSE)
    # } else if (modellist_g$cdist == "gaussian") {
    #     fy <- dnorm(Y_data, muy, sqrt(epsilon))
    # } else if (modellist_g$cdist == "nbinom") {
    #     fy <- dnbinom(Y_data, mu = muy, size = epsilon)
    # } else if (modellist_g$cdist == "binomial") {
    #     fy <- dbinom(Y_data, size = 1, prob = muy)
    # }


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

    no_lv <- 0L

    # 2. Estimation of the Predictor Model
    if (!is.null(pmodel)) {
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
        lvnames <- vnames[vnames %in% lavaan::lavNames(fit_pm, "lv")]
        no_lv <- length(lvnames)
        # manifest (stochastic) variables
        ovnames <- vnames[vnames %in% lavaan::lavNames(fit_pm, "ov")]
        # manifest (fixed) variables
        fxnames <- vnames[!vnames %in% c(lvnames, ovnames)]
        # manifest variables (fixed & stochastic)
        opvnames <- c(ovnames, fxnames)
    } else {
        lvnames <- ovnames <- NULL
        fxnames <- opvnames <- vnames
        fit_pm <- NULL
    }

    if (no_lv) {
        # 4. Computation of the person-specific Gauss-Hermite points
        # for numerical integration
        # formula from Wikipedia
        # https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
        # Eaton, Morris L. (1983). Multivariate Statistics: a Vector Space Approach. John Wiley and Sons. pp. 116–117

        # as far as I understand, we use regression factor scores and their variance
        # easiest way should be:
        mubar <- lavaan::lavPredict(
            fit_pm,
            method = "regression",
            acov = TRUE,
            assemble = TRUE
        )
        att_tmp <- attr(mubar, "acov")
        mubar <- mubar[, c(lvnames)] |> as.data.frame()
        attr(mubar, "acov") <- lapply(
            att_tmp,
            function(v) {
                v <- v[lvnames, lvnames] |> as.matrix()
                colnames(v) <- rownames(v) <- lvnames
                return(v)
            }
        )

        # Initialize GH grid
        # Pruning not meaningful for one dimension
        # (would remove all but one value)

        if (no_lv > 1) {
            prune <- TRUE
        } else if (no_lv == 1) {
            prune <- FALSE
        }

        GHgrid <- lavacreg:::creg_init_grid(
            Q = length(lvnames),
            ip = Nip,
            type = inttype,
            prune = prune
        )

        GH_i <- data.frame()
        for (i in 1:N) {
            level_i <- which(levels(groupvar) %in% groupvar[i])
            mubar_i <- mubar[i, ] |> unlist()
            sigbar_i <- attr(mubar, "acov")[[level_i]]

            # Apply individual mean structure (unit-specific)
            GHgrid_i <- lavacreg:::creg_adapt_grid(mubar_i, sigbar_i, GHgrid)
            colnames(GHgrid_i$X) <- lvnames

            GH_i <- rbind(GH_i, cbind(id = i, GHgrid_i$X, GHw = GHgrid_i$W))
        }

        dataext <- cbind(data[GH_i$id, ], GH_i)
        groupvarext <- groupvar[GH_i$id]
    } else {
        dataext <- data
        dataext$id <- 1:N
        groupvarext <- groupvar
    }
    mf <- model.frame(formula, dataext)
    mt <- attr(mf, "terms") # allow model.frame to have updated it


    Y <- model.response(mf, "any") # e.g. factors are allowed
    ## avoid problems with 1D arrays, but keep names
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) names(Y) <- nm
    }
    ## null model support
    X <- if (!is.empty.model(mt)) {
        model.matrix(mt, mf, contrasts = NULL)
    } else {
        matrix(, NROW(Y), 0L)
    }


    lvglm_objective <- function(x) {
        if (disp_logical) {
            coef <- matrix(x, ncol = ngroups)
            disppar <- coef[ncoef_g, ] |> as.numeric()
            coef <- coef[-ncoef_g, ] |> as.matrix(ncol = ngroups)
        } else {
            coef <- matrix(x, ncol = ngroups)
            disppar <- numeric()
        }


        linkinv <- family$linkinv
        dfun <- family$dfun
        dy <- numeric(nrow(dataext))

        for (g in 1:ngroups) {
            sel_g_id <- groupvarext == levels(groupvar)[g]
            eta <- (X[sel_g_id, ] %*% coef[, g])
            mu <- linkinv(eta)
            dy[sel_g_id] <- dfun(Y[sel_g_id], mu, disppar)
        }

        dataext$GHw * dy
        data_t <- data.frame(
            id = dataext$id,
            dy = dy
        )
        tmp <- setDT(data_t)[, .(x = sum(dy)), by = id]

        if (anyNA(tmp)) {
            return(+Inf)
        }
        # 6.6 Logarithmize individual likelihoods and sum
        out <- tmp$x |>
            log() |>
            sum()

        # Scaling and inversion for nlmnib
        obj <- (-1 * out) / N
        return(obj)
    }

    #
    # 7. Set starting values
    disp_logical <- FALSE
    if (is.na(family$dispersion)) {
        disp_logical <- TRUE
        ncoef_g <- ncol(X) + 1
        ncoef <- ncoef_g * ngroups
        coefstart <- rep(c(rep(0, ncol(X)), 1), ngroups)
    } else {
        ncoef_g <- ncol(X)
        ncoef <- ncoef_g * ngroups
        coefstart <- rep(0, ncoef)
    }


    # 8. Optimization of likelihood
    fit_om <- nlminb(
        start = coefstart,
        objective = lvglm_objective
    )

    if (fit_om$converegence == 0) {}


    if (se == "none") {
        vcov_fit <- matrix()
    } else {
        # Compute standard errors
        information <- numDeriv::hessian(
            func = lvglm_objective,
            x = fit_om$par
        )
        vcov_fit <- solve(information) / N
    }


    input <- list(
        N_groups = ngroups,
        N_cells = n_cell,
        N_coef_g = ncoef_g,
        N_coef = ncoef,
        prednames = colnames(X),
        formula = formula,
        se = se,
        family = family
    )


    # Return some results
    structure(
        list(
            call = cl,
            input = input,
            modelmatrix = X,
            dataext = dataext,
            groupvarext = groupvarext,
            groupvar = groupvar,
            fit_om = fit_om,
            fit_pm = fit_pm,
            vcov = vcov_fit
        ),
        class = "lvglm"
    )
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
        linkinv = stats$linkinv, dispersion = NA
    ), class = "family")
}
