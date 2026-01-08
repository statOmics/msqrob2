data(cptac)
cptac <- zeroIsNA(cptac, "peptides")
cptac <- logTransform(cptac, "peptides", "peptides_log")
cptac <- normalize(cptac, "peptides_log", "peptides_norm", "center.median")

test_that("msqrobLmer", {
    y <- assay(cptac[["peptides_norm"]])
    coldata <- colData(cptac)
    rowdata <- rowData(cptac)[["peptides"]]
    ## Default usage
    formula <- ~ condition + (1|lab)
    out <- msqrobLmer(y, formula, coldata)
    expect_identical(length(out), nrow(y))
    expect_true(all(sapply(out, class) == "StatModel"))
    expect_identical(
        as.list(table(sapply(out, function(x) x@type))),
        list(fitError = 104L, lmer = 252L)
    )
    expect_equal(
        out[[1]]@params,
        list(
            coefficients = c(
                `(Intercept)` = -1.43239555985396, conditionB = 0.303816995791296, conditionC = 0.400933593097784, conditionD = -0.716021732465039, `(Intercept)lablab2` = 0.277916037800659, `(Intercept)lablab3` = -0.277916037800659),
            vcovUnscaled = new(
                "dsCMatrix",
                i = c(0L, 0L, 1L, 0L, 1L, 2L, 0L, 1L, 2L, 3L, 0L, 1L, 2L, 3L, 4L, 0L, 1L, 2L, 3L, 4L, 5L),
                p = c(0L, 1L, 3L, 6L, 10L, 15L, 21L),
                Dim = c(6L, 6L),
                Dimnames = list(
                    c("(Intercept)", "conditionB", "conditionC", "conditionD", "(Intercept)lablab2", "(Intercept)lablab3"),
                    c("(Intercept)", "conditionB", "conditionC", "conditionD", "(Intercept)lablab2", "(Intercept)lablab3")
                ),
                x = c(0.972657631232327, -0.668226165745811, 1.16822616574581, -0.668226165745811, 0.668226165745811, 1.28592310873781, -0.668226165745811, 0.668226165745811, 0.651459442390789, 1.24184307090335, -0.304431465486516, 6.48495793590282e-17, -0.117696942992002, 0.0167667233550216, 0.422128408478519, -0.304431465486516, -1.3687424862664e-16, 0.117696942992002, -0.0167667233550217, 0.186734522494514, 0.422128408478519),
                uplo = "U", factors = list()
            ),
            sigma = 0.643133756500329,
            df.residual = 4.38661228982676,
            w = c(0.748249658021019, 0.748249658021019, 1, 1, 1, 1, 0.353602245028265, 0.397011079589108, 1)
        )
    )
    expect_equal( ## test spike-in coefs
        out[[grep("ups", rowdata$Proteins)[[1]]]]@params$coefficients[["conditionB"]],
        log2(0.74) - log2(0.25),
        tolerance = 0.3
    )
    expect_equal( ## test constant coefs
        out[[grep("YEAST", rowdata$Proteins)[[1]]]]@params$coefficients[["conditionB"]],
        0,
        tolerance = 0.5
    )
    ## Test with rowData annotations
    formula <- ~ condition + (1|lab) + (1 | Sequence)
    featureGroups <- rowdata$Proteins
    out <- msqrobLmer(
        y = y, formula = formula, data = coldata, rowdata = rowdata,
        featureGroups = featureGroups
    )
    expect_identical(
        as.list(table(sapply(out, function(x) x@type))),
        list(fitError = 14L, lmer = 46L)
    )
    expect_equal( ## test spike-in coefs
        out[[grep("ups", names(out))[[1]]]]@params$coefficients[["conditionB"]],
        log2(0.74) - log2(0.25),
        tolerance = 0.3
    )
    expect_equal( ## test constant coefs
        out[[grep("YEAST", names(out))[[1]]]]@params$coefficients[["conditionB"]],
        0,
        tolerance = 0.3
    )
    ## Disable robust
    out <- msqrobLmer(
        y = y, formula = formula, data = coldata, rowdata = rowdata,
        featureGroups = featureGroups, robust = FALSE
    )
    expect_identical(
        as.list(table(sapply(out, function(x) x@type))),
        list(fitError = 14L, lmer = 46L)
    )
    expect_equal( ## test spike-in coefs
        out[[grep("ups", names(out))[[1]]]]@params$coefficients[["conditionB"]],
        log2(0.74) - log2(0.25),
        tolerance = 0.3
    )
    expect_equal( ## test constant coefs
        out[[grep("YEAST", names(out))[[1]]]]@params$coefficients[["conditionB"]],
        0,
        tolerance = 0.3
    )
    expect_true(all(sapply(out, function(x) all(x@params$w == 1))))
    ## Enable ridge
    out <- msqrobLmer(
        y = y, formula = formula, data = coldata, rowdata = rowdata,
        featureGroups = featureGroups, ridge = TRUE
    )
    expect_identical(
        as.list(table(sapply(out, function(x) x@type))),
        list(fitError = 14L, lmer = 46L)
    )
    expect_equal( ## test spike-in coefs
        out[[grep("ups", names(out))[[1]]]]@params$coefficients[["ridgeconditionB"]],
        log2(0.74) - log2(0.25),
        tolerance = 0.3
    )
    expect_equal( ## test constant coefs
        out[[grep("YEAST", names(out))[[1]]]]@params$coefficients[["ridgeconditionB"]],
        0,
        tolerance = 0.1 ## decrease tolerance for non DA with ridge
    )
    ## Disable QR
    ## no impact if no ridge
    expect_identical(
        msqrobLmer(
            y = y, formula = formula, data = coldata, rowdata = rowdata,
            featureGroups = featureGroups
        ),
        msqrobLmer(
            y = y, formula = formula, data = coldata, rowdata = rowdata,
            featureGroups = featureGroups, doQR = FALSE
        )
    )
    ## disable QR when ridge
    out <- msqrobLmer(
        y = y, formula = formula, data = coldata, rowdata = rowdata,
        featureGroups = featureGroups, ridge = TRUE, doQR = FALSE
    )
    expect_identical(
        as.list(table(sapply(out, function(x) x@type))),
        list(fitError = 14L, lmer = 46L)
    )
    expect_equal( ## test spike-in coefs
        out[[grep("ups", names(out))[[1]]]]@params$coefficients[["ridgeconditionB"]],
        log2(0.74) - log2(0.25),
        tolerance = 0.3
    )
    expect_equal( ## test constant coefs
        out[[grep("YEAST", names(out))[[1]]]]@params$coefficients[["ridgeconditionB"]],
        0,
        tolerance = 0.1 ## decrease tolerance for non DA with ridge
    )
    ## Test M-estimation parameters
    out <- msqrobLmer(
        y = y, formula = formula, data = coldata, rowdata = rowdata,
        featureGroups = featureGroups, maxitRob = 10, tol = 1E-2
    )
    expect_identical(
        as.list(table(sapply(out, function(x) x@type))),
        list(fitError = 14L, lmer = 46L)
    )
    expect_equal( ## test spike-in coefs
        out[[grep("ups", names(out))[[1]]]]@params$coefficients[["conditionB"]],
        log2(0.74) - log2(0.25),
        tolerance = 0.3
    )
    expect_equal( ## test constant coefs
        out[[grep("YEAST", names(out))[[1]]]]@params$coefficients[["conditionB"]],
        0,
        tolerance = 0.3
    )
    ## Test lmer parameters
    stop("lmerArgs are not used...")
})
