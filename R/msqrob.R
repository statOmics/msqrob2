#' Function to fit msqrob models using lm and rlm
#'
#' @description Low-level function for parameter estimation with msqrob
#'              using the ordinary least squares or robust regression
#'              base on the MASS::rlm function.
#'
#' @param y A `matrix` with the quantified feature intensities. The
#'        features are along the rows and samples along the columns.
#'
#' @param formula Model formula. The model is built based on the
#'        covariates in the data object.
#'
#' @param data A `DataFrame` with information on the design. It has
#'        the same number of rows as the number of columns (samples) of
#'        `y`.
#'
#' @param robust `boolean(1)` to indicate if robust regression is
#'        performed to account for outliers. Default is `TRUE`. If
#'        `FALSE` an OLS fit is performed.
#'
#' @param maxitRob `numeric(1)` indicating the maximum iterations in
#'        the IRWLS algorithm used in the M-estimation step of the robust
#'        regression.
#'
#' @param trend `character(1)` or `FALSE` controlling the empirical Bayes
#'        variance trend. `FALSE` (default) or `"none"` uses a global prior
#'        variance. `"mean"` conditions the prior on mean log-intensity (like
#'        \code{limma::eBayes(trend = TRUE)}). `"count"` conditions it on
#'        log2 precursor count (like DEqMS). `"combined"` uses a linear
#'        projection of `log(s^2)` on both mean and log2 count as a 1-D
#'        covariate for \code{squeezeVar}.
#'
#' @param counts An optional numeric vector of length \code{nrow(y)} with the
#'        number of precursors (PSMs or peptides) per feature. Only used when
#'        \code{trend = "count"} or \code{trend = "combined"}.
#'
#' @examples
#'
#' # Load example data
#' # The data are a Feature object with containing
#' # a SummarizedExperiment named "peptide" with MaxQuant peptide intensities
#' # The data are a subset of spike-in the human-ecoli study
#' # The variable condition in the colData of the Feature object
#' # contains information on the spike in condition a-e (from low to high)
#' data(pe)
#'
#' # Aggregate peptide intensities in protein expression values
#' pe <- aggregateFeatures(pe, i = "peptide", fcol = "Proteins", name = "protein")
#' pe
#'
#' # Fit MSqrob model using robust regression with the MASS rlm function
#' models <- msqrobLm(assay(pe[["protein"]]), ~condition, colData(pe))
#' #' getCoef(models[[1]])
#' @return A list of objects of the `StatModel` class.
#'
#' @rdname msqrobLm
#'
#' @author Lieven Clement, Oliver M. Crook
#'
#' @importFrom MASS rlm
#' @importFrom stats model.matrix lm.fit
#' @importFrom limma squeezeVar
#' @importFrom methods is
#' @importFrom BiocParallel bplapply
#' @importFrom estimability nonest.basis all.estble
#'
#' @export
msqrobLm <- function(y,
    formula,
    data,
    robust = TRUE,
    maxitRob = 5,
    trend = FALSE,
    counts = NULL) {
    if (isFALSE(trend)) trend <- "none"
    trend <- match.arg(trend, c("none", "mean", "count", "combined"))
    myDesign <- model.matrix(formula, data)
    models <- BiocParallel::bplapply(asplit(y, 1),
        function(y, design) {
            obs <- is.finite(y)
            type <- "fitError"
            model <- list(
                coefficients = NA, vcovUnscaled = NA,
                sigma = NA, df.residual = NA, w = NA
            )
            
            # nb is the non-estimable basis calculated with the estimability 
            # package used to catch and solve issues with changes of the 
            # reference class due to missingness. 
            nb <- NULL
            if (sum(obs) > 0) {
                td <- .trim_design(design[obs, , drop = FALSE])
                X  <- td$X
                nb <- td$nb
                colnames_orig <- td$colnames_orig
                y <- y[obs]

                if (robust) {
                    mod <- try(MASS::rlm(X, y, method = "M", maxit = maxitRob),
                        silent = TRUE
                    )
                    if (!is(mod, "try-error")) type <- "rlm"
                } else {
                    mod <- try(lm.fit(X, y))
                    if (!is(mod, "try-error") && mod$rank == ncol(X)) type <- "lm"
                }

                if (type == "rlm") {
                    w <- mod$w
                    sw <- sum(w)
                    df.residual <- sw - mod$rank
                    sigma <- sqrt(sum(w * mod$resid^2) / df.residual)
                    if (df.residual < 2L) type <- "fitError"
                }

                if (type == "lm") {
                    w <- NULL
                    sigma <- sqrt(sum(mod$residuals^2) / mod$df.residual)
                    df.residual <- mod$df.residual
                    if (df.residual < 2L) type <- "fitError"
                }

                if (type != "fitError") {
                    coef <- rep(NA, length(colnames_orig))
                    names(coef) <- colnames_orig
                    coef[names(mod$coef)] <- mod$coef

                    vcovUnscaled <- matrix(NA,
                        nrow = length(colnames_orig),
                        ncol = length(colnames_orig)
                    )
                    rownames(vcovUnscaled) <- colnames(vcovUnscaled) <- colnames_orig
                    vcovUnscaled[names(mod$coef), names(mod$coef)] <- .vcovUnscaled(mod)

                    model <- list(
                        coefficients = coef,
                        vcovUnscaled = vcovUnscaled,
                        sigma = sigma,
                        df.residual = df.residual,
                        w = w
                    )
                }
            }
            StatModel(
                type = type,
                params = c(model, list(robust = robust, nonest.basis = nb)),
                varPosterior = as.numeric(NA),
                dfPosterior = as.numeric(NA)
            )
        },
        design = myDesign
    )

    covariate <- NULL
    if (trend != "none") {
        vars      <- vapply(models, getVar, numeric(1))
        mean_expr <- rowMeans(y, na.rm = TRUE)
        covariate <- switch(trend,
            mean     = .make_trend_covariate(vars, mean_expr = mean_expr),
            count    = if (!is.null(counts)) {
                .make_trend_covariate(vars, counts = counts)
            } else {
                warning("trend='count' requires counts; no trend applied.")
                NULL
            },
            combined = if (!is.null(counts)) {
                .make_trend_covariate(vars, mean_expr = mean_expr, counts = counts)
            } else {
                warning("trend='combined' without counts; falling back to 'mean'.")
                .make_trend_covariate(vars, mean_expr = mean_expr)
            }
        )
    }

    .squeezeAndUpdatePosteriors(models, covariate = covariate)
}



#' Function to fit msqrob models with ridge regression and/or random effects using lme4
#'
#' @description Low-level function for parameter estimation with msqrob
#'              using the robust ridge regression. The models can be fitted for each
#'              feature (e.g. summarised protein expression values) or multiple features
#'              belonging to the same accession can be modelled simultaneously
#'              e.g. peptide-based models where all peptide intensities for the same
#'              protein are modelled simultaneously. The fold changes and uncertainty
#'              estimates are then calculated at the protein level while correcting
#'              for peptide species and within sample correlation.
#'
#' @param y A `matrix` with the quantified feature intensities. The
#'        features are along the rows and samples along the columns.
#'
#' @param formula Model formula. The model is built based on the
#'        covariates in the data object.
#'
#' @param data A `DataFrame` with information on the design. It has
#'        the same number of rows as the number of columns (samples) of
#'        `y`.
#'
#' @param rowdata A `DataFrame` with the rowData information of the SummarizedExperiment.
#'        It has the same number of rows as the number of rows (features) of
#'        `y`.
#'
#' @param robust `boolean(1)` to indicate if robust regression is
#'        performed to account for outliers. Default is `TRUE`. If
#'        `FALSE` an OLS fit is performed.
#'
#' @param ridge `boolean(1)` to indicate if ridge regression is
#'        performed. Default is `FALSE`. If `TRUE` the fixed effects are
#'        estimated via penalized regression and shrunken to zero.
#'
#' @param maxitRob `numeric(1)` indicating the maximum iterations in
#'        the IRWLS algorithm used in the M-estimation step of the robust
#'        regression.
#'
#' @param tol `numeric(1)` indicating the tolerance for declaring convergence
#'        of the M-estimation loop.
#'
#'
#' @param doQR `boolean(1)` to indicate if a QR decomposition is applied to the
#'        fixed-effect design matrix before the Scheipl ridge encoding. Default is `TRUE`.
#'        When `TRUE` the predictor columns are orthogonalised so that shrinkage is
#'        invariant to predictor ordering and collinearity, and fixed effects are shrunken
#'        toward zero in the space of treatment contrasts. When `FALSE` the raw design
#'        matrix columns are used directly (standard L2 penalty on original parameters).
#'
#' @param featureGroups vector of type `character` or vector of type `factor` indicating how to aggregate
#'        the features. Is only used when multiple features are used to build the model, e.g. when starting
#'        from peptide data and modelling the fold change at the protein level. The default is `NULL`
#'
#' @param trend `character(1)` or `FALSE` controlling the empirical Bayes
#'        variance trend. `FALSE` (default) or `"none"` uses a global prior
#'        variance. `"mean"` conditions the prior on mean log-intensity. `"count"`
#'        conditions it on log2 precursor count. `"combined"` uses a linear
#'        projection of `log(s^2)` on both as a 1-D covariate for
#'        \code{squeezeVar}. When \code{featureGroups} groups multiple features
#'        per protein, precursor counts are inferred automatically; for
#'        protein-level data supply them via \code{counts}.
#'
#' @param counts An optional numeric vector (length = number of unique feature
#'        groups) giving the number of precursors per protein. Used when
#'        \code{trend = "count"} or \code{trend = "combined"} and each feature
#'        group contains a single row (i.e. the data are already at protein level).
#'
#' @param lmerArgs a list (of correct class, resulting from 'lmerControl()'
#'        containing control parameters, including the nonlinear optimizer to be used
#'        and parameters to be passed through to the nonlinear optimizer, see the
#'        'lmerControl' documentation of the lme4 package for more details.
#'        Default is `list(control = lmerControl(calc.derivs = FALSE))`
#'
#' @examples
#'
#' # Load example data
#' # The data are a Feature object with containing
#' # a SummarizedExperiment named "peptide" with MaxQuant peptide intensities
#' # The data are a subset of spike-in the human-ecoli study
#' # The variable condition in the colData of the Feature object
#' # contains information on the spike in condition a-e (from low to high)
#' data(pe)
#'
#' # Aggregate peptide intensities in protein expression values
#' pe <- aggregateFeatures(pe, i = "peptide", fcol = "Proteins", name = "protein")
#'
#' # Fit MSqrob model using robust ridge regression upon summarization of
#' # peptide intensities into protein expression values
#' modelsRidge <- msqrobLmer(assay(pe[["protein"]]), ~condition, data = colData(pe),
#'                           ridge = TRUE)
#' getCoef(modelsRidge[[1]])
#'
#' # Fit MSqrob model using robust ridge regression starting from peptide intensities
#' # The fold changes are calculated at the protein level while correcting for
#' # the different peptide species in each sample and the correlation between
#' # peptide intensities of peptides of the same protein in the same sample.
#' # Add the samples variable to colData
#' colData(pe)$samples <- rownames(colData(pe))
#' modelsPepBased <- msqrobLmer(assay(pe[["peptide"]]),
#'     formula = ~condition + (1|samples) + (1|Sequence), data = colData(pe),
#'     rowdata = rowData(pe[["peptide"]]), featureGroups = rowData(pe[["peptide"]])$Proteins,
#'     ridge = TRUE)
#' getCoef(modelsPepBased[[1]])
#' @return A list of objects of the `StatModel` class.
#'
#' @rdname msqrobLmer
#'
#' @author Lieven Clement, Oliver M. Crook
#'
#' @importFrom MASS psi.huber
#' @importFrom stats resid update.formula resid mad
#' @importFrom methods as is
#' @import lme4
#' @import Matrix
#' @importFrom BiocParallel bplapply bpmapply
#' @importFrom MultiAssayExperiment DataFrame
#'
#' @export

msqrobLmer <- function(y,
    formula,
    data,
    rowdata = NULL,
    tol = 1e-6,
    robust = TRUE,
    ridge = FALSE,
    maxitRob = 1,
    doQR = TRUE,
    featureGroups = NULL,
    trend = FALSE,
    counts = NULL,
    lmerArgs = list(control = lmerControl(calc.derivs = FALSE))) {

    if (isFALSE(trend)) trend <- "none"
    trend <- match.arg(trend, c("none", "mean", "count", "combined"))

    if (is.null(featureGroups)) featureGroups <- rownames(y)

    if (!is.null(rowdata)) {
        rowdata <- rowdata[colnames(rowdata) %in% all.vars(formula)]
        rowdata <- split.data.frame(rowdata, featureGroups)
    }

    data <- data[, colnames(data) %in% all.vars(formula), drop = FALSE]
    y <- split.data.frame(y, featureGroups)

    ## Resolve the count vector for trend covariates that need it.
    ## Priority: user-supplied counts > inferred from featureGroups > none.
    if (trend %in% c("count", "combined")) {
        if (!is.null(counts)) {
            if (!is.null(names(counts))) counts <- counts[names(y)]
        } else {
            pep_counts <- sapply(y, nrow)
            if (all(pep_counts == 1L)) {
                warning(
                    "trend='", trend, "' requires precursor counts but all ",
                    "feature groups contain a single feature. ",
                    "Supply counts= for protein-level data. ",
                    if (trend == "combined") "Falling back to 'mean'."
                    else "No trend applied."
                )
                counts <- NULL
                trend  <- if (trend == "combined") "mean" else "none"
            } else {
                counts <- pep_counts
            }
        }
    }

    mean_expr <- if (trend %in% c("mean", "combined")) {
        sapply(y, function(yi) mean(as.matrix(yi), na.rm = TRUE))
    } else NULL

    FUN <- if (ridge) .ridge_msqrobLmer else .noridge_msqrobLmer
    extraArgs <- if (ridge) {
        list(formula = formula, coldata = data, doQR = doQR,
             robust = robust, maxitRob = maxitRob, tol = tol,
             lmerArgs = lmerArgs)
    } else {
        list(formula = formula, coldata = data,
             robust = robust, maxitRob = maxitRob, tol = tol,
             lmerArgs = lmerArgs)
    }
    models <- .lmer_apply(y, rowdata, FUN, extraArgs)

    covariate <- if (trend != "none") {
        vars <- vapply(models, getVar, numeric(1))
        switch(trend,
            mean     = .make_trend_covariate(vars, mean_expr = mean_expr),
            count    = .make_trend_covariate(vars, counts = counts),
            combined = .make_trend_covariate(vars, mean_expr = mean_expr,
                                             counts = counts)
        )
    } else NULL

    .squeezeAndUpdatePosteriors(models, covariate = covariate)
}


#' Function to fit msqrob models to peptide counts using glm
#'
#' @description Low-level function for parameter estimation with msqrob
#'              by modeling peptide counts using quasibinomial glm
#'
#' @param y A `matrix` with the peptide counts. The
#'        features are along the rows and samples along the columns.
#'
#' @param npep A vector with number of peptides per protein. It has as length
#'        the number of rows of y. The counts are equal or larger than the largest
#'        peptide count in y.
#'
#' @param formula Model formula. The model is built based on the
#'        covariates in the data object.
#'
#' @param data A `DataFrame` with information on the design. It has
#'        the same number of rows as the number of columns (samples) of
#'        `y`.
#'
#' @param priorCount A 'numeric(1)', which is a prior count to be added to the observations to shrink
#'          the estimated log-fold-changes towards zero.
#'
#' @param binomialBound logical, if 'TRUE' then the quasibinomial variance estimator will
#'        be never smaller than 1 (no underdispersion).
#'
#' @examples
#'
#' # Load example data
#' # The data are a Feature object with containing
#' # a SummarizedExperiment named "peptide" with MaxQuant peptide intensities
#' # The data are a subset of spike-in the human-ecoli study
#' # The variable condition in the colData of the Feature object
#' # contains information on the spike in condition a-e (from low to high)
#' data(pe)
#'
#' # Aggregate peptide intensities in protein expression values
#' pe <- aggregateFeatures(pe, i = "peptide", fcol = "Proteins", name = "protein")
#' pe
#'
#' # Fit MSqrob model using robust regression with the MASS rlm function
#' models <- msqrobGlm(
#'     aggcounts(pe[["protein"]]),
#'     rowData(pe[["protein"]])[[".n"]],
#'     ~condition,
#'     colData(pe)
#' )
#' getCoef(models[[1]])
#' @return A list of objects of the `StatModel` class.
#'
#' @rdname msqrobGlm
#'
#' @author Lieven Clement
#'
#' @importFrom limma squeezeVar
#' @importFrom stats model.matrix glm.fit binomial
#' @importFrom methods is
#' @importFrom BiocParallel bplapply
#'
#' @export

msqrobGlm <- function(y,
    npep,
    formula,
    data,
    priorCount = .1,
    binomialBound = TRUE) {
    myDesign <- model.matrix(formula, data)
    models <- BiocParallel::bplapply(seq_len(nrow(y)),
        function(i, y, npep, myDesign) {
            type <- "fitError"
            model <- list(
                coefficients = NA, vcovUnscaled = NA,
                sigma = NA, df.residual = NA, w = NULL
            )
            mod <- NULL
            if (npep[i] >= max(y[i, ])) {
                mod <- try(glm.fit(
                    y = cbind(y[i, ], npep[i] - y[i, ]) + priorCount,
                    x = myDesign,
                    family = binomial()
                ))
                if (!is(mod, "try-error") && mod$rank == ncol(myDesign)) {
                    type <- "quasibinomial"
                }
            }
            if (!is.null(mod) && !is(mod, "try-error")) {
                if (mod$deviance < 0) {
                    mod$deviance <- sum(pmax(
                        mod$family$dev.resids(mod$y, mod$fitted.values, mod$prior.weights), 0
                    ))
                }
                if (mod$df.residual < 2L) type <- "fitError"
            }
            if (type != "fitError") {
                model <- list(
                    coefficients = mod$coef,
                    vcovUnscaled = .vcovUnscaled(mod),
                    sigma = sqrt(mod$deviance / mod$df.residual),
                    df.residual = mod$df.residual,
                    w = mod$w
                )
            }
            StatModel(
                type = type,
                params = model,
                varPosterior = as.numeric(NA),
                dfPosterior = as.numeric(NA)
            )
        },
        y = y, npep = npep, myDesign = myDesign
    )

    .squeezeAndUpdatePosteriors(models, binomialBound = binomialBound)
}


## Build a 1-D trend covariate for squeezeVar.
## mean_expr and/or counts control which predictors are included:
##   one predictor  → pass it directly (squeezeVar's loess handles the rest)
##   both           → linear projection of log(s²) ~ mean + log2(count) onto a
##                    1-D axis on the same scale squeezeVar uses internally.
## The QR rank-reduction handles constant count vectors (e.g. all features have
## one precursor) by silently dropping the redundant column.
.make_trend_covariate <- function(vars, mean_expr = NULL, counts = NULL) {
    if (is.null(mean_expr) && is.null(counts)) return(NULL)
    if (is.null(counts))    return(mean_expr)
    if (is.null(mean_expr)) return(log2(pmax(counts, 0.5)))
    log_s2 <- ifelse(is.finite(vars) & vars > 0, log(vars), NA_real_)
    ok <- is.finite(log_s2) #& is.finite(mean_expr) & is.finite(log2(pmax(counts, 0.5)))
    # If we have less than 3 variances that are finite we cannot fit a model with
    # an intercept and two slope parameters
    if (sum(ok) < 3L) return(mean_expr)
    xmat <- cbind(1, mean_expr, log2(pmax(counts, 0.5)))
    #qrX  <- qr(xmat, tol = 1e-10)
    #if (qrX$rank < ncol(xmat))
    #    xmat <- xmat[, qrX$pivot[seq_len(qrX$rank)], drop = FALSE]
    fit           <- lm.fit(xmat[ok,], log_s2[ok])
    #covariate     <- rep(NA_real_, length(vars))
    covariate <- xmat %*% fit$coef
    return(covariate)
}


## Squeeze sample variances together via empirical Bayes posterior means and
## update varPosterior / dfPosterior slots in-place.
.squeezeAndUpdatePosteriors <- function(models, binomialBound = FALSE,
                                        covariate = NULL) {
    hlp <- limma::squeezeVar(
        var       = vapply(models, getVar, numeric(1)),
        df        = vapply(models, getDF,  numeric(1)),
        covariate = covariate
    )
    for (i in seq_along(models)) {
        vp <- as.numeric(hlp$var.post[i])
        df <- as.numeric(hlp$df.prior + getDF(models[[i]]))
        if (binomialBound && !is.na(vp) && vp < 1) {
            vp <- 1
            df <- Inf
        }
        models[[i]]@varPosterior <- vp
        models[[i]]@dfPosterior  <- df
    }
    models
}


## Dispatch bplapply or bpmapply depending on whether rowdata is provided.
.lmer_apply <- function(y, rowdata, FUN, args) {
    if (is.null(rowdata)) {
        do.call(BiocParallel::bplapply, c(list(X = y, FUN = FUN), args))
    } else {
        do.call(BiocParallel::bpmapply, c(list(FUN = FUN, y, rowdata), list(MoreArgs = args)))
    }
}


## Drop all-zero columns, compute nonest.basis, and rank-reduce the design matrix.
## Returns X (trimmed), colnames_orig, nb, and rank_deficient (TRUE only when QR
## reduction removed columns beyond the all-zero pass).
.trim_design <- function(X) {
    colnames_orig <- colnames(X)
    nb <- estimability::nonest.basis(X)
    if (!identical(nb, estimability::all.estble)) rownames(nb) <- colnames_orig
    X <- X[, colMeans(X == 0) != 1, drop = FALSE]
    qrX <- qr(X)
    rank_deficient <- qrX$rank < ncol(X)
    if (rank_deficient)
        X <- X[, qrX$pivot[seq_len(qrX$rank)], drop = FALSE]
    list(X = X, colnames_orig = colnames_orig, nb = nb,
         rank_deficient = rank_deficient)
}


## Pad coefficients and vcovUnscaled back to the full colnames_orig space with NAs.
.pad_to_orig_space <- function(params, colnames_orig) {
    ran_names <- names(params$coefficients)[
        !names(params$coefficients) %in% colnames_orig]
    all_names <- c(colnames_orig, ran_names)
    coef_full <- rep(NA_real_, length(all_names))
    names(coef_full) <- all_names
    coef_full[names(params$coefficients)] <- params$coefficients
    vcov_full <- matrix(NA_real_, length(all_names), length(all_names))
    rownames(vcov_full) <- colnames(vcov_full) <- all_names
    vcov_full[names(params$coefficients), names(params$coefficients)] <-
        params$vcovUnscaled
    params$coefficients <- coef_full
    params$vcovUnscaled <- vcov_full
    params
}


## Initialise weights, run optional robust IRWLS, and extract model components.
.extract_lmer_fit <- function(model, robust, maxitRob, tol) {
    model@frame$`(weights)` <- rep(1, nrow(model@frame))
    sseOld <- model@devcomp$cmp["pwrss"]
    if (robust) model <- .robust_fitting(model, maxitRob, sseOld, tol)
    df <- .getDfLmer(model)
    list(
        model        = model,
        sigma        = sigma(model),
        betas        = .getBetaB(model),
        vcovUnscaled = as.matrix(.getVcovBetaBUnscaled(model)),
        df.residual  = if (is.na(df)) 0 else df
    )
}


## Fit the mixed models with ridge regression
.ridge_msqrobLmer <- function(y, rowdata = NULL, formula, coldata,
    doQR = TRUE, robust, maxitRob = 1, tol = 1e-06,
    lmerArgs = list(control = lmerControl(calc.derivs = FALSE))) {

    data <- .create_data(y, rowdata, coldata)

    fixed  <- model.matrix(nobars(formula), data = data)
    data$y <- as.matrix(y)
    obs    <- as.vector(!is.na(data$y))
    data   <- data[obs, , drop = FALSE]
    td     <- .trim_design(fixed[obs, , drop = FALSE])
    colnames_orig <- td$colnames_orig
    nb            <- td$nb
    data$fixed    <- td$X

    has_intercept <- colnames(data$fixed)[1] == "(Intercept)"

    if (sum(!grepl("(Intercept)", colnames(fixed))) < 2 &&
            !identical(nobars(formula)[[2]], 1)) {
        stop(
            "The mean model must have more than two parameters for ridge regression.\n",
            "If you really want to adopt ridge regression when your factor has only two levels\n",
            "rerun the function with a formula where you drop the intercept. e.g. ~-1+condition"
        )
    }

    ## Build the Z matrix for the Scheipl ridge random effect
    if (doQR) {
        qrFixed <- qr(data$fixed)
        Zridge  <- qr.Q(qrFixed)
        colnames(Zridge) <- colnames(data$fixed)
        if (has_intercept) Zridge <- Zridge[, -1, drop = FALSE]
    } else {
        Zridge <- if (has_intercept) data$fixed[, -1, drop = FALSE] else data$fixed
    }

    if (is.null(findbars(formula))) {
        formula <- formula(y ~ (1 | ridge))
    } else {
        if (!identical(nobars(formula)[[2]], 1)) {
            formula <- formula(paste0(
                "y ~ (1|ridge) + ",
                paste0("(", paste(findbars(formula), collapse = ")+("), ")")
            ))
        } else {
            formula <- update.formula(formula, y ~ .)
        }
    }

    model <- NULL
    try({
        data$ridge <- factor(rep(colnames(Zridge), length = nrow(data)), levels = colnames(Zridge))

        parsedFormulaC <- lFormula(formula, data = as.list(data))
        parsedFormulaC$reTrms$cnms$ridge <- ""
        ridgeId <- grep(names(parsedFormulaC$reTrms$Ztlist), pattern = "ridge")
        parsedFormulaC$reTrms$Ztlist[[ridgeId]] <-
            as(Matrix(t(Zridge)), class(parsedFormulaC$reTrms$Ztlist[[ridgeId]]))
        parsedFormulaC$reTrms$Zt <- do.call(rbind, parsedFormulaC$reTrms$Ztlist)

        devianceFunctionC <- do.call(mkLmerDevfun, parsedFormulaC)
        optimizerOutputC  <- optimizeLmer(devianceFunctionC)
        model <- mkMerMod(
            rho    = environment(devianceFunctionC),
            opt    = optimizerOutputC,
            reTrms = parsedFormulaC$reTrms,
            fr     = parsedFormulaC$fr
        )
    }, silent = TRUE)

    type   <- "fitError"
    params <- list(coefficients = NA, vcovUnscaled = NA, sigma = NA, df.residual = NA, w = NA)

    if (!is.null(model)) {
        try({
            fit          <- .extract_lmer_fit(model, robust, maxitRob, tol)
            betas        <- fit$betas
            vcovUnscaled <- fit$vcovUnscaled
            coefNames    <- names(betas)

            ## Back-transform from QR space to original parameter space
            if (doQR && any(grepl("ridge", names(betas)))) {
                ids <- if (has_intercept) c(1L, grep("ridge", names(betas))) else grep("ridge", names(betas))
                Rinv <- diag(length(betas))
                Rinv[ids, ids] <- solve(qr.R(qrFixed))
                if (has_intercept) Rinv[1L, 1L] <- 1
                betas        <- c(Rinv %*% betas)
                names(betas) <- coefNames
                vcovUnscaled <- Rinv %*% vcovUnscaled %*% t(Rinv)
                rownames(vcovUnscaled) <- colnames(vcovUnscaled) <- names(betas)
            }

            ## Rename ridge BLUPs from internal "ridgeconditionX" names to the
            ## original design-matrix parameter names so getContrast can find them.
            ridge_ids <- grep("ridge", names(betas))
            if (length(ridge_ids) > 0L) {
                names(betas)[ridge_ids]           <- colnames(Zridge)
                rownames(vcovUnscaled)[ridge_ids] <- colnames(Zridge)
                colnames(vcovUnscaled)[ridge_ids] <- colnames(Zridge)
            }

            params <- .create_model(betas, vcovUnscaled, fit$sigma, fit$df.residual, fit$model)
            if (!is.na(params$df.residual)) {
                type   <- "lmer"
                params <- .pad_to_orig_space(params, colnames_orig)
            }
        }, silent = TRUE)
    }

    StatModel(
        type         = type,
        params       = c(params, list(robust = robust, ridge = TRUE, doQR = doQR,
                                      lmerArgs = lmerArgs, nonest.basis = nb)),
        varPosterior = as.numeric(NA),
        dfPosterior  = as.numeric(NA)
    )
}


## Fit the mixed models without ridge regression
.noridge_msqrobLmer <- function(y, rowdata = NULL, formula, coldata,
    robust, maxitRob = 0, tol = 1e-06,
    lmerArgs = list(control = lmerControl(calc.derivs = FALSE))) {

    data <- .create_data(y, rowdata, coldata)

    mm <- model.matrix(nobars(formula), data = data)
    formula <- update.formula(formula, y ~ .)

    data$y <- as.matrix(y)
    obs    <- as.vector(!is.na(data$y))
    data   <- data[obs, , drop = FALSE]
    td     <- .trim_design(mm[obs, , drop = FALSE])
    colnames_orig     <- td$colnames_orig
    nb                <- td$nb
    data_model_matrix <- td$X

    if (td$rank_deficient) {
        for (nm in colnames(data_model_matrix))
            data[[paste0(".x.", nm)]] <- data_model_matrix[, nm]
        fix     <- paste(paste0("`.x.", colnames(data_model_matrix), "`"), collapse = " + ")
        bars    <- findbars(formula)
        formula <- if (is.null(bars)) {
            as.formula(paste("y ~", fix))
        } else {
            as.formula(paste("y ~", fix, "+",
                             paste0("(", paste(bars, collapse = ")+("), ")")))
        }
    }

    model <- NULL
    try({
        model <- lmer(formula, as.data.frame(data))
    }, silent = TRUE)

    type   <- "fitError"
    params <- list(coefficients = NA, vcovUnscaled = NA, sigma = NA, df.residual = NA, w = NA)

    if (!is.null(model)) {
        try({
            fit    <- .extract_lmer_fit(model, robust, maxitRob, tol)
            params <- .create_model(fit$betas, fit$vcovUnscaled, fit$sigma, fit$df.residual, fit$model)
            if (!is.na(params$df.residual)) {
                type   <- "lmer"
                params <- .pad_to_orig_space(params, colnames_orig)
            }
        }, silent = TRUE)
    }

    StatModel(
        type         = type,
        params       = c(params, list(robust = robust, ridge = FALSE, lmerArgs = lmerArgs,
                                      nonest.basis = nb)),
        varPosterior = as.numeric(NA),
        dfPosterior  = as.numeric(NA)
    )
}


## Calculate unscaled covariance matrix for lm or rlm fit
.vcovUnscaled <- function(model) {
    p1  <- seq_len(model$rank)
    p   <- length(model$coefficients)
    out <- matrix(NA, p, p)
    out[!is.na(model$coefficients), !is.na(model$coefficients)] <-
        chol2inv(model$qr$qr[p1, p1, drop = FALSE])
    colnames(out) <- rownames(out) <- names(model$coefficients)
    out
}

#' @import purrr
.zNames <- function(model) {
    ranefLevels <- purrr::imap(model@flist, ~ paste0(.y, levels(.x)))
    unlist(lapply(seq_along(model@cnms),
        function(x, cnms, levels) c(outer(cnms[[x]], levels[[names(cnms)[x]]], paste0)),
        cnms = model@cnms, levels = ranefLevels
    ))
}

#' @import lme4
#' @import Matrix
#' @importFrom methods cbind2

.getVcovBetaBUnscaled <- function(model) {
    X  <- lme4::getME(model, "X")
    Z  <- lme4::getME(model, "Z")
    XZ <- cbind2(X, Z)

    if (is.null(model@frame$`(weights)`)) model@frame$`(weights)` <- 1

    vcovInv <- Matrix::crossprod(model@frame$`(weights)`^.5 * XZ)
    Ginv <- Matrix::solve(
        Matrix::tcrossprod(getME(model, "Lambda")) +
            Matrix::Diagonal(ncol(Z), 1e-18)
    )

    i <- -seq_len(ncol(X))
    vcovInv[i, i] <- vcovInv[i, i] + Ginv
    vcovInv <- Matrix::solve(vcovInv)

    rownames(vcovInv) <- colnames(vcovInv) <- c(colnames(X), .zNames(model))
    vcovInv
}

#' @import lme4
.getBetaB <- function(model) {
    betaB <- c(
        as.vector(lme4::getME(model, "beta")),
        as.vector(lme4::getME(model, "b"))
    )
    names(betaB) <- c(colnames(model@pp$X), .zNames(model))
    betaB
}

#' @importFrom stats resid
.getDfLmer <- function(object) {
    w <- object@frame$"(weights)"
    if (is.null(w)) w <- 1
    sigma <- sigma(object)
    sum((resid(object) * sqrt(w))^2) / sigma^2
}

.robust_fitting <- function(model, maxitRob, sseOld, tol) {
    while (maxitRob > 0) {
        maxitRob <- maxitRob - 1
        res <- resid(model)
        model@frame$`(weights)` <- MASS::psi.huber(res / mad(res, 0))
        model <- refit(model)
        sse <- model@devcomp$cmp["pwrss"]
        if (abs(sseOld - sse) / sseOld <= tol) break
        sseOld <- sse
    }
    model
}

.create_model <- function(betas, vcovUnscaled, sigma, df.residual, model) {
    if (df.residual < 2L) {
        list(coefficients = NA, vcovUnscaled = NA, sigma = NA, df.residual = NA, w = NA)
    } else {
        list(
            coefficients = betas,
            vcovUnscaled = vcovUnscaled,
            sigma        = sigma,
            df.residual  = df.residual,
            w            = model@frame$`(weights)`
        )
    }
}

.create_data <- function(y, rowdata, coldata) {
    nr <- nrow(coldata)
    if (is.null(rowdata)) {
        coldata[rep(seq_len(nr), each = nrow(y)), , drop = FALSE]
    } else {
        data <- cbind(
            coldata[rep(seq_len(nr), each = nrow(y)), ],
            rowdata[rep(seq_len(nrow(rowdata)), ncol(y)), ]
        )
        data <- DataFrame(data)
        colnames(data) <- c(colnames(coldata), colnames(rowdata))
        data
    }
}
