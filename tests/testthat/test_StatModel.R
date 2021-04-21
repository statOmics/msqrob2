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
