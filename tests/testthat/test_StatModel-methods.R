.create_minimal_data <- function(){
    reference_present_no_ref <- c(FALSE, FALSE)
    names(reference_present_no_ref) <- c("conditionb", "conditionc")
    
    reference_present_ref <- c(TRUE, TRUE)
    names(reference_present_ref) <- c("conditionb", "conditionc")

    vcovu_no_ref <- matrix(c(0.2, -0.2, NA, -0.2, 0.4, NA, NA, NA, NA), nrow = 3, byrow = TRUE)
    colnames(vcovu_no_ref) <- c("(Intercept)", "conditionb", "conditionc")
    rownames(vcovu_no_ref) <- c("(Intercept)", "conditionb", "conditionc")

    vcovu_ref <- matrix(c(0.2, -0.2, -0.2, -0.2, 0.4, -0.2, -0.2, -0.2, 0.6 ), nrow = 3, byrow = TRUE)
    colnames(vcovu_ref) <- c("(Intercept)", "conditionb", "conditionc")
    rownames(vcovu_ref) <- c("(Intercept)", "conditionb", "conditionc")

    stat_model_list <- list(feat1 = StatModel(
        type = "lm",
        params = list(coefficients = c("(Intercept)" = 1, "conditionb" = 2, "conditionc" = NA),
            df.residual = 10, sigma = 0.5,
            vcovUnscaled = vcovu_no_ref,
            referencePresent = reference_present_no_ref),
        varPosterior = 0.1,
        dfPosterior = 6
    ), feat2 = StatModel(
        type = "lm",
        params = list(coefficients = c("(Intercept)" = 1, "conditionb" = 2, "conditionc" = 3),
            df.residual = 10, sigma = 0.5,
            vcovUnscaled = vcovu_ref,
            referencePresent = reference_present_ref),
        varPosterior = 0.1,
        dfPosterior = 6
    ))
    return(stat_model_list)
}

test_that("getContrast", {
    stat_model_list <- .create_minimal_data()
    L <- matrix(c(1, 0), nrow = 2, byrow = TRUE)
    colnames(L) <- "conditionb=0"
    rownames(L) <- c("conditionb", "conditionc")

    res_NA <- matrix(NA)
    rownames(res_NA) <- "conditionb=0"
    res_NA_dbl <- matrix(as.double(NA))
    rownames(res_NA_dbl) <- "conditionb=0"
    expect_equal(res_NA, getContrast(stat_model_list$feat1, L))
    expect_equal(res_NA_dbl, getContrast(stat_model_list$feat1, L, TRUE))

    res_2 <- matrix(2)
    rownames(res_2) <- "conditionb=0"

    L_no_0 <- L[L != 0, , drop = FALSE] # What occurs in topFeatures
    expect_identical(res_2, getContrast(stat_model_list$feat1, L_no_0, TRUE))

    expect_identical(res_2, getContrast(stat_model_list$feat2, L))
    expect_identical(res_2, getContrast(stat_model_list$feat2, L, TRUE))
})

test_that("varContrast", {
    stat_model_list <- .create_minimal_data()
    L <- matrix(c(1, 0), nrow = 2, byrow = TRUE)
    colnames(L) <- "conditionb=0"
    rownames(L) <- c("conditionb", "conditionc")

    res_NA <- matrix(NA)
    rownames(res_NA) <- "conditionb=0"
    colnames(res_NA) <- "conditionb=0"

    res_NA_dbl <- matrix(as.double(NA))
    rownames(res_NA_dbl) <- "conditionb=0"
    colnames(res_NA_dbl) <- "conditionb=0"
    expect_equal(res_NA, varContrast(stat_model_list$feat1, L))
    expect_equal(res_NA_dbl, varContrast(stat_model_list$feat1, L, TRUE))

    res_2 <- matrix(0.04)
    rownames(res_2) <- "conditionb=0"
    colnames(res_2) <- "conditionb=0"
    L_no_0 <- L[L != 0, , drop = FALSE] # What occurs in topFeatures
    expect_equal(res_2, varContrast(stat_model_list$feat1, L_no_0, TRUE))

    expect_equal(res_2, varContrast(stat_model_list$feat2, L))
    expect_equal(res_2, varContrast(stat_model_list$feat2, L, TRUE))
})
