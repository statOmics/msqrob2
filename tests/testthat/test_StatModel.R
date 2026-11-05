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
