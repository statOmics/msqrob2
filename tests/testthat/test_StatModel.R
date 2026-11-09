test_that("StatModel construction", {
    ## Empty mode
    expect_true(validObject(StatModel()))
    ## Default model type
    expect_identical(StatModel()@type, "fitError")
    ## Fully defined model
    mod <- StatModel(
        type = "rlm",
        params = list(x = 3, y = 7, b = 4),
        varPosterior = c(0.1, 0.2, 0.3),
        dfPosterior = c(6, 7, 8)
    )
    expect_true(validObject(mod))
})

test_that("getContrast and varContrast return NA for failed fits", {
    ## a failed fit as msqrobLm / msqrobLmer store it
    mod <- StatModel(
        type = "fitError",
        params = list(coefficients = NA, vcovUnscaled = NA, sigma = NA,
                      df.residual = NA, w = NA),
        varPosterior = as.numeric(NA),
        dfPosterior = as.numeric(NA)
    )
    L <- makeContrast("conditionB = 0", parameterNames = "conditionB")
    expect_true(is.na(getContrast(mod, L)[1, 1]))
    expect_true(is.na(varContrast(mod, L)[1, 1]))
})

test_that("getContrast and varContrast return NA for names that are no model parameters", {
    set.seed(1)
    cd <- data.frame(condition = factor(rep(c("A", "B"), each = 4)))
    y  <- matrix(rnorm(8), 1, dimnames = list("p1", NULL))
    mod <- msqrobLm(y, ~condition, cd, robust = FALSE)[[1]]
    ## "conditionX" is a typo: only the first contrast can be computed
    L <- makeContrast(c("conditionB = 0", "conditionX = 0", "conditionB - conditionX = 0"),
                      parameterNames = c("conditionB", "conditionX"))
    est <- getContrast(mod, L)
    expect_equal(unname(est[1, 1]), getCoef(mod)[["conditionB"]])
    expect_true(all(is.na(est[2:3, 1])))
    v <- varContrast(mod, L)
    expect_equal(unname(v[1, 1]), getVcovUnscaled(mod)["conditionB", "conditionB"] * getVarPosterior(mod))
    expect_true(all(is.na(v[2:3, ])))
    expect_true(all(is.na(v[, 2:3])))
})
