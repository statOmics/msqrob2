## Tests for msqrob.R: msqrobLm, msqrobGlm, msqrobLmer, and internal helpers

library(testthat)
library(msqrob2)
library(BiocParallel)
library(S4Vectors)
library(MASS)
library(lme4)

register(SerialParam())

# ---------------------------------------------------------------------------
# Shared toy-data builders
# ---------------------------------------------------------------------------

.make_lm_data <- function(n_feat = 5, n_rep = 4, seed = 42) {
    set.seed(seed)
    n_samples <- 2L * n_rep
    cond <- factor(rep(c("A", "B"), each = n_rep))
    coldata <- DataFrame(condition = cond,
                         row.names = paste0("s", seq_len(n_samples)))
    y <- matrix(rnorm(n_feat * n_samples, mean = 10, sd = 1),
                nrow = n_feat, ncol = n_samples,
                dimnames = list(paste0("feat", seq_len(n_feat)),
                                rownames(coldata)))
    # Add a small effect in condition B for all but the last feature
    y[-n_feat, cond == "B"] <- y[-n_feat, cond == "B"] + 2
    list(y = y, coldata = coldata)
}

## Repeated-measures design: each subject measured in both conditions.
## This makes ~condition + (1|subject) identifiable for lmer.
.make_lmer_data <- function(n_feat = 4, n_subj = 5, seed = 7) {
    set.seed(seed)
    n_samples <- 2L * n_subj   # each subject appears twice
    cond    <- factor(rep(c("A", "B"), times = n_subj))
    subject <- factor(rep(paste0("subj", seq_len(n_subj)), each = 2))
    snames  <- paste0("s", seq_len(n_samples))
    coldata <- DataFrame(condition = cond, subject = subject,
                         row.names = snames)
    y <- matrix(rnorm(n_feat * n_samples, mean = 10, sd = 1),
                nrow = n_feat, ncol = n_samples,
                dimnames = list(paste0("feat", seq_len(n_feat)), snames))
    y[, cond == "B"] <- y[, cond == "B"] + 2
    list(y = y, coldata = coldata)
}

## 3-condition design required for ridge regression (>2 mean-model parameters).
.make_ridge_data <- function(n_feat = 4, n_rep = 4, seed = 99) {
    set.seed(seed)
    n_cond    <- 3L
    n_samples <- n_cond * n_rep
    cond    <- factor(rep(paste0("cond", seq_len(n_cond)), each = n_rep))
    coldata <- DataFrame(condition = cond,
                         row.names = paste0("s", seq_len(n_samples)))
    y <- matrix(rnorm(n_feat * n_samples, mean = 10, sd = 1),
                nrow = n_feat, ncol = n_samples,
                dimnames = list(paste0("feat", seq_len(n_feat)),
                                rownames(coldata)))
    list(y = y, coldata = coldata)
}


# ===========================================================================
# msqrobLm
# ===========================================================================

test_that("msqrobLm returns a list of StatModel objects of correct length", {
    d <- .make_lm_data()
    models <- msqrobLm(d$y, ~condition, d$coldata)
    expect_type(models, "list")
    expect_length(models, nrow(d$y))
    expect_true(all(vapply(models, is, logical(1), "StatModel")))
})

test_that("msqrobLm robust=TRUE produces 'rlm' fit types for complete data", {
    d <- .make_lm_data()
    models <- msqrobLm(d$y, ~condition, d$coldata, robust = TRUE)
    types <- vapply(models, getFitMethod, character(1))
    expect_true(all(types %in% c("rlm", "fitError")))
    # With 4 replicates per group the majority should converge
    expect_gt(sum(types == "rlm"), 0L)
})

test_that("msqrobLm robust=FALSE produces 'lm' fit types for complete data", {
    d <- .make_lm_data()
    models <- msqrobLm(d$y, ~condition, d$coldata, robust = FALSE)
    types <- vapply(models, getFitMethod, character(1))
    expect_true(all(types %in% c("lm", "fitError")))
    expect_gt(sum(types == "lm"), 0L)
})

test_that("msqrobLm coefficients are named by design-matrix columns", {
    d <- .make_lm_data()
    models <- msqrobLm(d$y, ~condition, d$coldata, robust = FALSE)
    expected_names <- colnames(model.matrix(~condition, d$coldata))
    valid <- Filter(function(m) getFitMethod(m) != "fitError", models)
    expect_gt(length(valid), 0L)
    for (m in valid) {
        expect_named(getCoef(m), expected_names)
    }
})

test_that("msqrobLm: vcovUnscaled is square and symmetric for 'lm' fits", {
    d <- .make_lm_data()
    models <- msqrobLm(d$y, ~condition, d$coldata, robust = FALSE)
    valid <- Filter(function(m) getFitMethod(m) == "lm", models)
    expect_gt(length(valid), 0L)
    for (m in valid) {
        V <- getVcovUnscaled(m)
        expect_equal(nrow(V), ncol(V))
        expect_equal(V, t(V), tolerance = 1e-10)
    }
})

test_that("msqrobLm: varPosterior and dfPosterior are non-NA after squeezing", {
    d <- .make_lm_data()
    models <- msqrobLm(d$y, ~condition, d$coldata)
    valid <- Filter(function(m) getFitMethod(m) != "fitError", models)
    for (m in valid) {
        expect_false(is.na(getVarPosterior(m)))
        expect_false(is.na(getDfPosterior(m)))
        expect_gt(getVarPosterior(m), 0)
    }
})

test_that("msqrobLm: all-NA row yields fitError", {
    d <- .make_lm_data()
    d$y[1, ] <- NA_real_
    models <- msqrobLm(d$y, ~condition, d$coldata)
    expect_equal(getFitMethod(models[[1]]), "fitError")
})

test_that("msqrobLm: row with too few observations yields fitError", {
    d <- .make_lm_data()
    # Keep only 2 finite values so df.residual < 2 after fitting ~condition (2 params)
    d$y[2, ] <- NA_real_
    d$y[2, 1:2] <- c(10, 11)
    models <- msqrobLm(d$y, ~condition, d$coldata, robust = FALSE)
    expect_equal(getFitMethod(models[[2]]), "fitError")
})

test_that("msqrobLm sigma is positive for valid fits", {
    d <- .make_lm_data()
    models <- msqrobLm(d$y, ~condition, d$coldata)
    valid <- Filter(function(m) getFitMethod(m) != "fitError", models)
    for (m in valid) {
        expect_gt(getSigma(m), 0)
    }
})

test_that("msqrobLm maxitRob argument is accepted without error", {
    d <- .make_lm_data()
    expect_no_error(msqrobLm(d$y, ~condition, d$coldata, robust = TRUE, maxitRob = 1))
    expect_no_error(msqrobLm(d$y, ~condition, d$coldata, robust = TRUE, maxitRob = 10))
})


# ===========================================================================
# msqrobGlm
# ===========================================================================

.make_glm_data <- function(n_prot = 5, n_rep = 4, seed = 13) {
    set.seed(seed)
    n_samples <- 2L * n_rep
    cond <- factor(rep(c("A", "B"), each = n_rep))
    coldata <- DataFrame(condition = cond,
                         row.names = paste0("s", seq_len(n_samples)))
    npep <- sample(3:8, n_prot, replace = TRUE)
    y <- t(vapply(npep, function(n) rbinom(n_samples, n, 0.6), integer(n_samples)))
    dimnames(y) <- list(paste0("prot", seq_len(n_prot)), rownames(coldata))
    list(y = y, npep = npep, coldata = coldata)
}

test_that("msqrobGlm returns a list of StatModel objects of correct length", {
    d <- .make_glm_data()
    models <- msqrobGlm(d$y, d$npep, ~condition, d$coldata)
    expect_type(models, "list")
    expect_length(models, nrow(d$y))
    expect_true(all(vapply(models, is, logical(1), "StatModel")))
})

test_that("msqrobGlm produces 'quasibinomial' type for valid counts", {
    d <- .make_glm_data()
    models <- msqrobGlm(d$y, d$npep, ~condition, d$coldata)
    types <- vapply(models, getFitMethod, character(1))
    expect_true(all(types %in% c("quasibinomial", "fitError")))
    expect_gt(sum(types == "quasibinomial"), 0L)
})

test_that("msqrobGlm: row with counts exceeding npep yields fitError", {
    d <- .make_glm_data()
    # Force a violation: set counts above npep for row 1
    d$y[1, ] <- d$npep[1] + 1L
    models <- msqrobGlm(d$y, d$npep, ~condition, d$coldata)
    expect_equal(getFitMethod(models[[1]]), "fitError")
})

test_that("msqrobGlm: varPosterior and dfPosterior are updated", {
    d <- .make_glm_data()
    models <- msqrobGlm(d$y, d$npep, ~condition, d$coldata)
    valid <- Filter(function(m) getFitMethod(m) != "fitError", models)
    for (m in valid) {
        expect_false(is.na(getVarPosterior(m)))
    }
})

test_that("msqrobGlm binomialBound=TRUE enforces varPosterior >= 1", {
    d <- .make_glm_data()
    models <- msqrobGlm(d$y, d$npep, ~condition, d$coldata, binomialBound = TRUE)
    valid <- Filter(function(m) getFitMethod(m) != "fitError", models)
    for (m in valid) {
        vp <- getVarPosterior(m)
        if (!is.na(vp)) expect_gte(vp, 1)
    }
})

test_that("msqrobGlm coefficients are named by design-matrix columns", {
    d <- .make_glm_data()
    models <- msqrobGlm(d$y, d$npep, ~condition, d$coldata)
    expected_names <- colnames(model.matrix(~condition, d$coldata))
    valid <- Filter(function(m) getFitMethod(m) == "quasibinomial", models)
    for (m in valid) {
        expect_named(getCoef(m), expected_names)
    }
})


# ===========================================================================
# msqrobLmer
# ===========================================================================

test_that("msqrobLmer ridge=FALSE returns list of StatModel objects", {
    d <- .make_lmer_data()
    models <- msqrobLmer(d$y, ~condition + (1 | subject),
                         data = d$coldata, ridge = FALSE)
    expect_type(models, "list")
    expect_length(models, nrow(d$y))
    expect_true(all(vapply(models, is, logical(1), "StatModel")))
})

test_that("msqrobLmer ridge=FALSE produces 'lmer' types for complete data", {
    d <- .make_lmer_data()
    models <- msqrobLmer(d$y, ~condition + (1 | subject),
                         data = d$coldata, ridge = FALSE)
    types <- vapply(models, getFitMethod, character(1))
    expect_true(all(types %in% c("lmer", "fitError")))
    expect_gt(sum(types == "lmer"), 0L)
})

test_that("msqrobLmer ridge=FALSE: varPosterior and dfPosterior updated", {
    d <- .make_lmer_data()
    models <- msqrobLmer(d$y, ~condition + (1 | subject),
                         data = d$coldata, ridge = FALSE)
    valid <- Filter(function(m) getFitMethod(m) == "lmer", models)
    expect_gt(length(valid), 0L)
    for (m in valid) {
        expect_false(is.na(getVarPosterior(m)))
        expect_false(is.na(getDfPosterior(m)))
    }
})

test_that("msqrobLmer ridge=TRUE doQR=TRUE runs without error", {
    d <- .make_ridge_data()
    expect_no_error(
        msqrobLmer(d$y, ~condition, data = d$coldata,
                   ridge = TRUE, doQR = TRUE)
    )
})

test_that("msqrobLmer ridge=TRUE doQR=FALSE runs without error", {
    d <- .make_ridge_data()
    expect_no_error(
        msqrobLmer(d$y, ~condition, data = d$coldata,
                   ridge = TRUE, doQR = FALSE)
    )
})

test_that("msqrobLmer ridge=TRUE produces 'lmer' types for most features", {
    d <- .make_ridge_data()
    models <- msqrobLmer(d$y, ~condition, data = d$coldata, ridge = TRUE)
    types <- vapply(models, getFitMethod, character(1))
    expect_true(all(types %in% c("lmer", "fitError")))
    expect_gt(sum(types == "lmer"), 0L)
})

test_that("msqrobLmer ridge=TRUE doQR=TRUE and doQR=FALSE fit same number of features", {
    d <- .make_ridge_data()
    mQR   <- msqrobLmer(d$y, ~condition, data = d$coldata, ridge = TRUE, doQR = TRUE)
    mNone <- msqrobLmer(d$y, ~condition, data = d$coldata, ridge = TRUE, doQR = FALSE)
    typesQR   <- vapply(mQR,   getFitMethod, character(1))
    typesNone <- vapply(mNone, getFitMethod, character(1))
    expect_equal(sum(typesQR == "lmer"), sum(typesNone == "lmer"))
})

test_that("msqrobLmer varPosterior and dfPosterior are updated for valid lmer fits", {
    d <- .make_ridge_data()
    models <- msqrobLmer(d$y, ~condition, data = d$coldata, ridge = TRUE)
    valid <- Filter(function(m) getFitMethod(m) == "lmer", models)
    expect_gt(length(valid), 0L)
    for (m in valid) {
        expect_false(is.na(getVarPosterior(m)))
        expect_false(is.na(getDfPosterior(m)))
    }
})

test_that("msqrobLmer ridge=TRUE with 2-level factor errors with informative message", {
    # 2-level factor with intercept is prohibited for ridge
    d <- .make_lm_data(n_feat = 3, n_rep = 3)
    expect_error(
        msqrobLmer(d$y, ~condition, data = d$coldata, ridge = TRUE),
        regexp = "more than two parameters"
    )
})

test_that("msqrobLmer featureGroups aggregates multiple peptides per protein", {
    d <- .make_lmer_data(n_feat = 6, n_subj = 4)
    groups <- factor(rep(c("prot1", "prot2"), each = 3))
    expect_no_error(
        models <- msqrobLmer(
            d$y,
            ~condition + (1 | subject),
            data          = d$coldata,
            featureGroups = groups,
            ridge         = FALSE
        )
    )
    expect_length(models, nlevels(groups))
})


# ===========================================================================
# Internal helpers (accessed via :::)
# ===========================================================================

## .vcovUnscaled -----------------------------------------------------------

test_that(".vcovUnscaled returns symmetric named matrix for lm.fit object", {
    set.seed(1)
    X <- cbind(1, rnorm(10), rnorm(10))
    colnames(X) <- c("(Intercept)", "x1", "x2")
    y <- X %*% c(5, 2, -1) + rnorm(10, sd = 0.5)
    fit <- lm.fit(X, y)
    V <- msqrob2:::.vcovUnscaled(fit)
    expect_equal(dim(V), c(ncol(X), ncol(X)))
    expect_equal(V, t(V), tolerance = 1e-10)
    expect_named(rownames(V), NULL)  # rownames present but not a named object
    expect_equal(colnames(V), colnames(X))
})

test_that(".vcovUnscaled returns NA-padded matrix for rank-deficient lm.fit", {
    set.seed(2)
    # Rank-deficient: x3 = x2
    X <- cbind(1, rnorm(10), rnorm(10))
    X <- cbind(X, X[, 2])
    colnames(X) <- c("(Intercept)", "x1", "x2", "x3")
    y <- X[, 1:3] %*% c(5, 2, -1) + rnorm(10, sd = 0.5)
    fit <- lm.fit(X, y)
    V <- msqrob2:::.vcovUnscaled(fit)
    # NA entries for the dropped column
    expect_true(any(is.na(V)))
    expect_equal(dim(V), c(ncol(X), ncol(X)))
})

## .create_data ------------------------------------------------------------

test_that(".create_data without rowdata stacks coldata nrow(y) times", {
    set.seed(3)
    nr <- 3L; ns <- 4L
    cd <- DataFrame(cond = factor(c("A","B","A","B")),
                    row.names = paste0("s", seq_len(ns)))
    y  <- matrix(rnorm(nr * ns), nr, ns)
    out <- msqrob2:::.create_data(y, NULL, cd)
    expect_equal(nrow(out), nr * ns)
    expect_equal(colnames(out), colnames(cd))
})

test_that(".create_data with rowdata adds rowdata columns", {
    set.seed(4)
    nr <- 2L; ns <- 3L
    cd <- DataFrame(cond = factor(c("A","B","A")),
                    row.names = paste0("s", seq_len(ns)))
    rd <- DataFrame(seq = c("PEP1", "PEP2"))
    y  <- matrix(rnorm(nr * ns), nr, ns)
    out <- msqrob2:::.create_data(y, rd, cd)
    expect_equal(nrow(out), nr * ns)
    expect_true("seq" %in% colnames(out))
})

## .create_model -----------------------------------------------------------

test_that(".create_model returns NA list when df.residual < 2", {
    dummy_model <- list(`(weights)` = rep(1, 3))
    # Create a minimal mock lmer frame
    fake_model <- list()
    fake_model$frame <- data.frame(`(weights)` = rep(1, 3))
    # Use a real lmer object to satisfy the interface
    # We test via the condition directly
    out <- msqrob2:::.create_model(
        betas        = c("(Intercept)" = 5, x = 2),
        vcovUnscaled = diag(2),
        sigma        = 1,
        df.residual  = 1,          # < 2 → should return NAs
        model        = list(frame = data.frame(`(weights)` = 1))
    )
    expect_true(is.na(out$coefficients))
    expect_true(is.na(out$sigma))
    expect_true(is.na(out$df.residual))
})

test_that(".create_model returns full list when df.residual >= 2", {
    # Need an actual lmer model for @frame access
    d <- .make_lmer_data()
    df <- as.data.frame(d$coldata)
    df$y <- as.numeric(d$y[1, ])
    lmer_mod <- lme4::lmer(y ~ condition + (1 | subject), data = df)
    betas <- c("(Intercept)" = 5, conditionB = 2)
    V <- diag(2)
    out <- msqrob2:::.create_model(
        betas        = betas,
        vcovUnscaled = V,
        sigma        = 1.2,
        df.residual  = 4,
        model        = lmer_mod
    )
    expect_equal(out$coefficients, betas)
    expect_equal(out$sigma, 1.2)
    expect_equal(out$df.residual, 4)
})

## .squeezeAndUpdatePosteriors ---------------------------------------------

test_that(".squeezeAndUpdatePosteriors updates varPosterior for all models", {
    # Build a list of StatModel objects with known sigma / df.residual
    mods <- lapply(seq_len(5), function(i) {
        StatModel(
            type   = "lm",
            params = list(coefficients = c("(Intercept)" = i),
                          vcovUnscaled = matrix(0.1),
                          sigma        = 0.5 + i * 0.1,
                          df.residual  = 5L,
                          w            = NULL)
        )
    })
    mods_squeezed <- msqrob2:::.squeezeAndUpdatePosteriors(mods)
    for (m in mods_squeezed) {
        expect_false(is.na(m@varPosterior))
        expect_false(is.na(m@dfPosterior))
        expect_gt(m@varPosterior, 0)
    }
})

test_that(".squeezeAndUpdatePosteriors with binomialBound enforces varPosterior >= 1", {
    # Force very small sigma so the squeezed variance would be < 1
    mods <- lapply(seq_len(5), function(i) {
        StatModel(
            type   = "quasibinomial",
            params = list(coefficients = c("(Intercept)" = 0.5),
                          vcovUnscaled = matrix(0.01),
                          sigma        = 0.001 * i,
                          df.residual  = 10L,
                          w            = NULL)
        )
    })
    mods_squeezed <- msqrob2:::.squeezeAndUpdatePosteriors(mods, binomialBound = TRUE)
    for (m in mods_squeezed) {
        vp <- m@varPosterior
        if (!is.na(vp)) expect_gte(vp, 1)
    }
})

test_that(".squeezeAndUpdatePosteriors preserves fitError models (NA sigma)", {
    good <- StatModel(
        type   = "lm",
        params = list(coefficients = 1, vcovUnscaled = matrix(0.1),
                      sigma = 1, df.residual = 4L, w = NULL)
    )
    bad <- StatModel(
        type   = "fitError",
        params = list(coefficients = NA, vcovUnscaled = NA,
                      sigma = NA, df.residual = NA, w = NA)
    )
    mods <- list(good, bad)
    out <- msqrob2:::.squeezeAndUpdatePosteriors(mods)
    # fitError model has NA varPosterior input; squeezeVar handles NaN/NA gracefully
    expect_true(is(out[[1]], "StatModel"))
    expect_true(is(out[[2]], "StatModel"))
})

## .lmer_apply -------------------------------------------------------------

test_that(".lmer_apply without rowdata calls FUN once per element of y", {
    counter <- 0L
    spy <- function(y, ...) { counter <<- counter + 1L; list() }
    y_list <- list(a = 1:3, b = 4:6, c = 7:9)
    msqrob2:::.lmer_apply(y_list, NULL, spy, list())
    expect_equal(counter, 3L)
})

test_that(".lmer_apply with rowdata calls FUN with both y and rowdata elements", {
    calls <- list()
    spy <- function(y, rd, ...) { calls[[length(calls) + 1L]] <<- list(y = y, rd = rd); list() }
    y_list  <- list(a = 1:3, b = 4:6)
    rd_list <- list(a = "rd_a", b = "rd_b")
    msqrob2:::.lmer_apply(y_list, rd_list, spy, list())
    expect_length(calls, 2L)
    expect_equal(calls[[1]]$rd, "rd_a")
    expect_equal(calls[[2]]$rd, "rd_b")
})


# ===========================================================================
# hypothesisTest: changing reference class
# ===========================================================================

## Build a 4-condition means-coded dataset.
## prot1: complete data.
## prot2: all condition-A replicates are NA (reference class absent).
## All six pairwise contrasts are tested; for prot2 the three comparisons
## that involve condition A must return NA, the other three must not.

.make_4cond_data <- function(n_rep = 4L, seed = 123L) {
    set.seed(seed)
    conds   <- factor(rep(c("A", "B", "C", "D"), each = n_rep))
    snames  <- paste0("s", seq_along(conds))
    coldata <- DataFrame(condition = conds, row.names = snames)
    y <- matrix(
        rnorm(2L * length(conds), mean = 10, sd = 1),
        nrow = 2L,
        dimnames = list(c("prot1", "prot2"), snames)
    )
    y["prot2", conds == "A"] <- NA_real_   # reference class absent for prot2
    list(y = y, coldata = coldata)
}

test_that("hypothesisTest: pairwise contrasts with absent reference class", {
    library(SummarizedExperiment)
    d  <- .make_4cond_data()
    se <- SummarizedExperiment(
        assays  = list(quant = d$y),
        colData = d$coldata
    )

    ## Fit using means coding so the absent condition is handled by the
    ## zero-column filter + NA-padding + estimability check.
    se <- msqrob(se, formula = ~0 + condition, overwrite = TRUE)

    ## All six pairwise contrasts among A, B, C, D.
    param <- c("conditionA", "conditionB", "conditionC", "conditionD")
    L <- makeContrast(
        c("conditionB - conditionA=0",
          "conditionC - conditionA=0",
          "conditionD - conditionA=0",
          "conditionB - conditionC=0",
          "conditionB - conditionD=0",
          "conditionC - conditionD=0"),
        parameterNames = param
    )
    se <- hypothesisTest(se, L)
    rd <- rowData(se)

    ## prot1 (complete data): all six contrasts estimable.
    for (contrast in colnames(L)) {
        logFC <- rd[[contrast]]["prot1", "logFC"]
        pval  <- rd[[contrast]]["prot1", "pval"]
        expect_false(is.na(logFC),
            info = paste("prot1 logFC should not be NA for contrast:", contrast))
        expect_false(is.na(pval),
            info = paste("prot1 pval should not be NA for contrast:", contrast))
    }

    ## prot2 (condA absent): contrasts involving A must be NA.
    ## Note: makeContrast strips the "=0" suffix from column names.
    a_contrasts     <- c("conditionB - conditionA",
                         "conditionC - conditionA",
                         "conditionD - conditionA")
    non_a_contrasts <- c("conditionB - conditionC",
                         "conditionB - conditionD",
                         "conditionC - conditionD")

    for (contrast in a_contrasts) {
        logFC <- rd[[contrast]]["prot2", "logFC"]
        pval  <- rd[[contrast]]["prot2", "pval"]
        expect_true(is.na(logFC),
            info = paste("prot2 logFC should be NA for contrast:", contrast))
        expect_true(is.na(pval),
            info = paste("prot2 pval should be NA for contrast:", contrast))
    }

    for (contrast in non_a_contrasts) {
        logFC <- rd[[contrast]]["prot2", "logFC"]
        pval  <- rd[[contrast]]["prot2", "pval"]
        expect_false(is.na(logFC),
            info = paste("prot2 logFC should not be NA for contrast:", contrast))
        expect_false(is.na(pval),
            info = paste("prot2 pval should not be NA for contrast:", contrast))
    }
})

test_that("hypothesisTest: pairwise contrasts with absent reference class, treatment coding", {
    library(SummarizedExperiment)
    d  <- .make_4cond_data()
    se <- SummarizedExperiment(
        assays  = list(quant = d$y),
        colData = d$coldata
    )

    ## Treatment coding: absent conditionA causes rank deficiency handled by
    ## pivot-drop + estimability check (not zero-column filter).
    se <- msqrob(se, formula = ~condition, overwrite = TRUE)

    ## Pairwise contrasts using treatment-coding parameter names.
    ## B-A, C-A, D-A are single coefficients; B-C, B-D, C-D are differences.
    param <- c("(Intercept)", "conditionB", "conditionC", "conditionD")
    L <- makeContrast(
        c("conditionB=0",
          "conditionC=0",
          "conditionD=0",
          "conditionB - conditionC=0",
          "conditionB - conditionD=0",
          "conditionC - conditionD=0"),
        parameterNames = param
    )
    se <- hypothesisTest(se, L)
    rd <- rowData(se)

    ## prot1 (complete data): all six contrasts estimable.
    for (contrast in colnames(L)) {
        logFC <- rd[[contrast]]["prot1", "logFC"]
        pval  <- rd[[contrast]]["prot1", "pval"]
        expect_false(is.na(logFC),
            info = paste("prot1 logFC should not be NA for contrast:", contrast))
        expect_false(is.na(pval),
            info = paste("prot1 pval should not be NA for contrast:", contrast))
    }

    ## prot2 (condA absent): A-involving contrasts non-estimable, others estimable.
    ## makeContrast strips "=0", so single-param contrasts keep the param name.
    a_contrasts     <- c("conditionB", "conditionC", "conditionD")
    non_a_contrasts <- c("conditionB - conditionC",
                         "conditionB - conditionD",
                         "conditionC - conditionD")

    for (contrast in a_contrasts) {
        logFC <- rd[[contrast]]["prot2", "logFC"]
        pval  <- rd[[contrast]]["prot2", "pval"]
        expect_true(is.na(logFC),
            info = paste("prot2 logFC should be NA for contrast:", contrast))
        expect_true(is.na(pval),
            info = paste("prot2 pval should be NA for contrast:", contrast))
    }

    for (contrast in non_a_contrasts) {
        logFC <- rd[[contrast]]["prot2", "logFC"]
        pval  <- rd[[contrast]]["prot2", "pval"]
        expect_false(is.na(logFC),
            info = paste("prot2 logFC should not be NA for contrast:", contrast))
        expect_false(is.na(pval),
            info = paste("prot2 pval should not be NA for contrast:", contrast))
    }
})

## ---------------------------------------------------------------------------
## Numerical accuracy: logFC == group mean difference for robust=FALSE
## ---------------------------------------------------------------------------

.expected_diffs_means <- function(y, cond) {
    mu <- sapply(levels(cond), function(lv)
              rowMeans(y[, cond == lv, drop = FALSE], na.rm = TRUE))
    list(
        "conditionB - conditionA" = mu[, "B"] - mu[, "A"],
        "conditionC - conditionA" = mu[, "C"] - mu[, "A"],
        "conditionD - conditionA" = mu[, "D"] - mu[, "A"],
        "conditionB - conditionC" = mu[, "B"] - mu[, "C"],
        "conditionB - conditionD" = mu[, "B"] - mu[, "D"],
        "conditionC - conditionD" = mu[, "C"] - mu[, "D"]
    )
}

.expected_diffs_treatment <- function(y, cond) {
    mu <- sapply(levels(cond), function(lv)
              rowMeans(y[, cond == lv, drop = FALSE], na.rm = TRUE))
    list(
        "conditionB"              = mu[, "B"] - mu[, "A"],
        "conditionC"              = mu[, "C"] - mu[, "A"],
        "conditionD"              = mu[, "D"] - mu[, "A"],
        "conditionB - conditionC" = mu[, "B"] - mu[, "C"],
        "conditionB - conditionD" = mu[, "B"] - mu[, "D"],
        "conditionC - conditionD" = mu[, "C"] - mu[, "D"]
    )
}

test_that("logFC equals group mean differences for means coding (robust=FALSE)", {
    library(SummarizedExperiment)
    d  <- .make_4cond_data()
    se <- SummarizedExperiment(
        assays  = list(quant = d$y),
        colData = d$coldata
    )
    se <- msqrob(se, formula = ~0 + condition, robust = FALSE, overwrite = TRUE)

    param <- c("conditionA", "conditionB", "conditionC", "conditionD")
    L <- makeContrast(
        c("conditionB - conditionA=0",
          "conditionC - conditionA=0",
          "conditionD - conditionA=0",
          "conditionB - conditionC=0",
          "conditionB - conditionD=0",
          "conditionC - conditionD=0"),
        parameterNames = param
    )
    se  <- hypothesisTest(se, L)
    rd  <- rowData(se)
    exp <- .expected_diffs_means(d$y, d$coldata$condition)

    ## prot1: all 6 contrasts estimable.
    for (nm in names(exp)) {
        expect_equal(rd[[nm]]["prot1", "logFC"], unname(exp[[nm]]["prot1"]),
            tolerance = 1e-10,
            info = paste("prot1 means-coding logFC mismatch:", nm))
    }

    ## prot2: conditionA absent — only non-A contrasts are estimable.
    non_a <- c("conditionB - conditionC", "conditionB - conditionD", "conditionC - conditionD")
    for (nm in non_a) {
        expect_equal(rd[[nm]]["prot2", "logFC"], unname(exp[[nm]]["prot2"]),
            tolerance = 1e-10,
            info = paste("prot2 means-coding logFC mismatch:", nm))
    }
})

test_that("logFC equals group mean differences for treatment coding (robust=FALSE)", {
    library(SummarizedExperiment)
    d  <- .make_4cond_data()
    se <- SummarizedExperiment(
        assays  = list(quant = d$y),
        colData = d$coldata
    )
    se <- msqrob(se, formula = ~condition, robust = FALSE, overwrite = TRUE)

    param <- c("(Intercept)", "conditionB", "conditionC", "conditionD")
    L <- makeContrast(
        c("conditionB=0",
          "conditionC=0",
          "conditionD=0",
          "conditionB - conditionC=0",
          "conditionB - conditionD=0",
          "conditionC - conditionD=0"),
        parameterNames = param
    )
    se  <- hypothesisTest(se, L)
    rd  <- rowData(se)
    exp <- .expected_diffs_treatment(d$y, d$coldata$condition)

    ## prot1: all 6 contrasts estimable.
    for (nm in names(exp)) {
        expect_equal(rd[[nm]]["prot1", "logFC"], unname(exp[[nm]]["prot1"]),
            tolerance = 1e-10,
            info = paste("prot1 treatment-coding logFC mismatch:", nm))
    }

    ## prot2: conditionA absent — only non-A contrasts are estimable.
    non_a <- c("conditionB - conditionC", "conditionB - conditionD", "conditionC - conditionD")
    for (nm in non_a) {
        expect_equal(rd[[nm]]["prot2", "logFC"], unname(exp[[nm]]["prot2"]),
            tolerance = 1e-10,
            info = paste("prot2 treatment-coding logFC mismatch:", nm))
    }
})

## ===========================================================================
## Heart design: location * tissue + patient (blocked factorial)
## ===========================================================================

## Synthetic 4-cell × 3-patient heart dataset with three missingness patterns.
## Cells: LA (left atrium), LV (left ventricle), RA (right atrium), RV (right ventricle)
## Contrasts of interest (matching the heart tutorial):
##   C1: tissueV                       = LV - LA effect (left region)
##   C2: tissueV + locationR:tissueV   = RV - RA effect (right region)
##   C3: tissueV + 0.5*locationR:tissueV = average V-A effect
##   C4: locationR:tissueV             = interaction (C2 - C1)
.make_heart_data <- function(seed = 42L) {
    set.seed(seed)
    loc <- factor(rep(c("L", "R"), each = 6), levels = c("L", "R"))
    tis <- factor(rep(rep(c("A", "V"), each = 3), 2), levels = c("A", "V"))
    pat <- factor(rep(c("3", "4", "8"), times = 4))
    sn  <- paste0(as.character(loc), as.character(tis), as.character(pat))
    cd  <- DataFrame(location = loc, tissue = tis, patient = pat, row.names = sn)
    y   <- matrix(
        rnorm(4L * 12L, mean = 10, sd = 1),
        nrow = 4L,
        dimnames = list(c("prot_cmp", "prot_LA", "prot_RV", "prot_p3"), sn)
    )
    y["prot_LA", cd$location == "L" & cd$tissue == "A"] <- NA_real_
    y["prot_RV", cd$location == "R" & cd$tissue == "V"] <- NA_real_
    y["prot_p3", cd$patient  == "3"]                    <- NA_real_
    list(y = y, coldata = cd)
}

## Expected logFC for the four heart contrasts when robust = FALSE.
## For balanced (or patient-missing) data, OLS estimates reduce to group mean
## differences; na.rm = TRUE averages over the observed patients only.
.heart_expected <- function(y, cd) {
    gm <- function(l, t)
        rowMeans(y[, cd$location == l & cd$tissue == t, drop = FALSE], na.rm = TRUE)
    C1 <- gm("L", "V") - gm("L", "A")
    C2 <- gm("R", "V") - gm("R", "A")
    list(C1 = C1, C2 = C2, C3 = 0.5 * (C1 + C2), C4 = C2 - C1)
}

test_that("heart design: logFC equals group mean contrasts (robust=FALSE)", {
    library(SummarizedExperiment)
    d  <- .make_heart_data()
    se <- SummarizedExperiment(assays = list(quant = d$y), colData = d$coldata)
    se <- msqrob(se, formula = ~ location * tissue + patient,
                 robust = FALSE, overwrite = TRUE)
    param <- colnames(model.matrix(~ location * tissue + patient,
                                   data = as.data.frame(d$coldata)))
    L <- makeContrast(
        c("tissueV = 0",
          "tissueV + locationR:tissueV = 0",
          "tissueV + 0.5*locationR:tissueV = 0",
          "locationR:tissueV = 0"),
        parameterNames = param
    )
    se  <- hypothesisTest(se, L)
    rd  <- rowData(se)
    cn  <- colnames(L)
    exp <- .heart_expected(d$y, d$coldata)

    ## prot_cmp: complete data, all 4 contrasts estimable.
    for (i in seq_along(cn)) {
        expect_equal(
            rd[[cn[i]]]["prot_cmp", "logFC"],
            unname(exp[[i]]["prot_cmp"]),
            tolerance = 1e-10,
            info = paste("prot_cmp heart contrast", i, "logFC mismatch:", cn[i])
        )
    }

    ## prot_p3: patient-3 absent, all 4 contrasts still estimable,
    ## logFC computed from patients 4 and 8 only.
    for (i in seq_along(cn)) {
        expect_equal(
            rd[[cn[i]]]["prot_p3", "logFC"],
            unname(exp[[i]]["prot_p3"]),
            tolerance = 1e-10,
            info = paste("prot_p3 heart contrast", i, "logFC mismatch:", cn[i])
        )
    }
})

test_that("heart design: estimability and logFC when a full cell is absent (robust=FALSE)", {
    library(SummarizedExperiment)
    d  <- .make_heart_data()
    se <- SummarizedExperiment(assays = list(quant = d$y), colData = d$coldata)
    se <- msqrob(se, formula = ~ location * tissue + patient,
                 robust = FALSE, overwrite = TRUE)
    param <- colnames(model.matrix(~ location * tissue + patient,
                                   data = as.data.frame(d$coldata)))
    L <- makeContrast(
        c("tissueV = 0",
          "tissueV + locationR:tissueV = 0",
          "tissueV + 0.5*locationR:tissueV = 0",
          "locationR:tissueV = 0"),
        parameterNames = param
    )
    se  <- hypothesisTest(se, L)
    rd  <- rowData(se)
    cn  <- colnames(L)
    exp <- .heart_expected(d$y, d$coldata)

    ## prot_LA (all LA absent): the intercept is non-estimable, so only
    ## C2 = RV - RA (which does not involve LA) remains estimable.
    expect_true (is.na(rd[[cn[1]]]["prot_LA", "logFC"]),
        info = "prot_LA C1 should be NA (LA absent)")
    expect_false(is.na(rd[[cn[2]]]["prot_LA", "logFC"]),
        info = "prot_LA C2 should not be NA (RV and RA both observed)")
    expect_true (is.na(rd[[cn[3]]]["prot_LA", "logFC"]),
        info = "prot_LA C3 should be NA (average involves non-estimable C1)")
    expect_true (is.na(rd[[cn[4]]]["prot_LA", "logFC"]),
        info = "prot_LA C4 should be NA (interaction involves non-estimable C1)")
    expect_equal(
        rd[[cn[2]]]["prot_LA", "logFC"],
        unname(exp$C2["prot_LA"]),
        tolerance = 1e-10,
        info = "prot_LA C2 logFC should equal mean(RV) - mean(RA)"
    )

    ## prot_RV (all RV absent): the interaction parameter is non-estimable,
    ## so only C1 = LV - LA (which does not involve RV) remains estimable.
    expect_false(is.na(rd[[cn[1]]]["prot_RV", "logFC"]),
        info = "prot_RV C1 should not be NA (LV and LA both observed)")
    expect_true (is.na(rd[[cn[2]]]["prot_RV", "logFC"]),
        info = "prot_RV C2 should be NA (RV absent)")
    expect_true (is.na(rd[[cn[3]]]["prot_RV", "logFC"]),
        info = "prot_RV C3 should be NA (average involves non-estimable C2)")
    expect_true (is.na(rd[[cn[4]]]["prot_RV", "logFC"]),
        info = "prot_RV C4 should be NA (interaction non-estimable without RV)")
    expect_equal(
        rd[[cn[1]]]["prot_RV", "logFC"],
        unname(exp$C1["prot_RV"]),
        tolerance = 1e-10,
        info = "prot_RV C1 logFC should equal mean(LV) - mean(LA)"
    )
})

## ===========================================================================
## rlm and ridge: NA pattern preserved under missingness
## ===========================================================================

test_that("4-cond means coding: NA pattern preserved for rlm and ridge", {
    library(SummarizedExperiment)
    d   <- .make_4cond_data()
    se  <- SummarizedExperiment(assays = list(quant = d$y), colData = d$coldata)
    param <- c("conditionA", "conditionB", "conditionC", "conditionD")
    L <- makeContrast(
        c("conditionB - conditionA=0",
          "conditionC - conditionA=0",
          "conditionD - conditionA=0",
          "conditionB - conditionC=0",
          "conditionB - conditionD=0",
          "conditionC - conditionD=0"),
        parameterNames = param
    )
    cn_A  <- colnames(L)[1:3]
    cn_nA <- colnames(L)[4:6]

    for (cfg in list(
        list(robust = TRUE,  ridge = FALSE, label = "rlm"),
        list(robust = TRUE,  ridge = TRUE,  label = "ridge")
    )) {
        se2 <- msqrob(se, formula = ~0 + condition,
                      robust = cfg$robust, ridge = cfg$ridge, overwrite = TRUE)
        se2 <- hypothesisTest(se2, L)
        rd  <- rowData(se2)
        lbl <- cfg$label

        ## prot1: all conditions present — all contrasts non-NA
        for (cn in colnames(L))
            expect_false(is.na(rd[[cn]]["prot1", "logFC"]),
                info = paste(lbl, "prot1", cn, "should be non-NA"))

        ## prot2: condA absent — A-involving contrasts NA, non-A contrasts non-NA
        for (cn in cn_A)
            expect_true(is.na(rd[[cn]]["prot2", "logFC"]),
                info = paste(lbl, "prot2", cn, "should be NA (condA absent)"))
        for (cn in cn_nA)
            expect_false(is.na(rd[[cn]]["prot2", "logFC"]),
                info = paste(lbl, "prot2", cn, "should be non-NA (condA not involved)"))
    }
})

test_that("4-cond treatment coding: NA pattern preserved for rlm and ridge", {
    library(SummarizedExperiment)
    d   <- .make_4cond_data()
    se  <- SummarizedExperiment(assays = list(quant = d$y), colData = d$coldata)
    param <- c("(Intercept)", "conditionB", "conditionC", "conditionD")
    L <- makeContrast(
        c("conditionB=0",
          "conditionC=0",
          "conditionD=0",
          "conditionB - conditionC=0",
          "conditionB - conditionD=0",
          "conditionC - conditionD=0"),
        parameterNames = param
    )
    cn_A  <- colnames(L)[1:3]
    cn_nA <- colnames(L)[4:6]

    for (cfg in list(
        list(robust = TRUE,  ridge = FALSE, label = "rlm"),
        list(robust = TRUE,  ridge = TRUE,  label = "ridge")
    )) {
        se2 <- msqrob(se, formula = ~condition,
                      robust = cfg$robust, ridge = cfg$ridge, overwrite = TRUE)
        se2 <- hypothesisTest(se2, L)
        rd  <- rowData(se2)
        lbl <- cfg$label

        ## prot1: all conditions present — all contrasts non-NA
        for (cn in colnames(L))
            expect_false(is.na(rd[[cn]]["prot1", "logFC"]),
                info = paste(lbl, "prot1", cn, "should be non-NA"))

        ## prot2: condA absent — A-involving contrasts NA, non-A contrasts non-NA
        for (cn in cn_A)
            expect_true(is.na(rd[[cn]]["prot2", "logFC"]),
                info = paste(lbl, "prot2", cn, "should be NA (condA absent)"))
        for (cn in cn_nA)
            expect_false(is.na(rd[[cn]]["prot2", "logFC"]),
                info = paste(lbl, "prot2", cn, "should be non-NA (condA not involved)"))
    }
})

test_that("heart design: estimability NA pattern preserved for rlm and ridge", {
    library(SummarizedExperiment)
    d   <- .make_heart_data()
    se  <- SummarizedExperiment(assays = list(quant = d$y), colData = d$coldata)
    param <- colnames(model.matrix(~ location * tissue + patient,
                                   data = as.data.frame(d$coldata)))
    L <- makeContrast(
        c("tissueV = 0",
          "tissueV + locationR:tissueV = 0",
          "tissueV + 0.5*locationR:tissueV = 0",
          "locationR:tissueV = 0"),
        parameterNames = param
    )
    cn <- colnames(L)

    for (cfg in list(
        list(robust = TRUE,  ridge = FALSE, label = "rlm"),
        list(robust = TRUE,  ridge = TRUE,  label = "ridge")
    )) {
        se2 <- msqrob(se, formula = ~ location * tissue + patient,
                      robust = cfg$robust, ridge = cfg$ridge, overwrite = TRUE)
        se2 <- hypothesisTest(se2, L)
        rd  <- rowData(se2)
        lbl <- cfg$label

        ## prot_cmp: all cells present, all contrasts non-NA
        for (i in seq_along(cn))
            expect_false(is.na(rd[[cn[i]]]["prot_cmp", "logFC"]),
                info = paste(lbl, "prot_cmp", cn[i], "should be non-NA"))

        ## prot_LA: LA absent — only C2 (RV - RA) estimable
        expect_true (is.na(rd[[cn[1]]]["prot_LA", "logFC"]),
            info = paste(lbl, "prot_LA C1 should be NA (LA absent)"))
        expect_false(is.na(rd[[cn[2]]]["prot_LA", "logFC"]),
            info = paste(lbl, "prot_LA C2 should be non-NA (RV and RA both present)"))
        expect_true (is.na(rd[[cn[3]]]["prot_LA", "logFC"]),
            info = paste(lbl, "prot_LA C3 should be NA"))
        expect_true (is.na(rd[[cn[4]]]["prot_LA", "logFC"]),
            info = paste(lbl, "prot_LA C4 should be NA"))

        ## prot_RV: RV absent — only C1 (LV - LA) estimable
        expect_false(is.na(rd[[cn[1]]]["prot_RV", "logFC"]),
            info = paste(lbl, "prot_RV C1 should be non-NA (LV and LA both present)"))
        expect_true (is.na(rd[[cn[2]]]["prot_RV", "logFC"]),
            info = paste(lbl, "prot_RV C2 should be NA (RV absent)"))
        expect_true (is.na(rd[[cn[3]]]["prot_RV", "logFC"]),
            info = paste(lbl, "prot_RV C3 should be NA"))
        expect_true (is.na(rd[[cn[4]]]["prot_RV", "logFC"]),
            info = paste(lbl, "prot_RV C4 should be NA"))

        ## prot_p3: patient-3 absent, all contrasts still estimable
        for (i in seq_along(cn))
            expect_false(is.na(rd[[cn[i]]]["prot_p3", "logFC"]),
                info = paste(lbl, "prot_p3", cn[i], "should be non-NA"))
    }
})
