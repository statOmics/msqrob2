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


test_that("topFeatures", {
    stat_model_list <- .create_minimal_data()
    L <- matrix(c(1, 0), nrow = 2, byrow = TRUE)
    colnames(L) <- "conditionb=0"
    rownames(L) <- c("conditionb", "conditionc")

    base_df <- data.frame("logFC" = c(2,2), "se" = c(0.2, 0.2),
        "df" = c(6,6), "t" = c(10,10),
        "pval" = c(5.791983e-05, 5.791983e-05),
        "adjPval" = c(5.791983e-05, 5.791983e-05),
        row.names = c("feat1", "feat2"))
    accept_df <- base_df
    accept_df$DifferentReference <- c(TRUE, FALSE)
    no_accept_df <- base_df[c(2,1),]
    no_accept_df["feat1", ] <- NA
    no_accept_df["feat1", "df"] <- 6

    expect_equal(no_accept_df,
        topFeatures(stat_model_list, L),
        tolerance = 1e-3)
    expect_equal(accept_df,
        topFeatures(stat_model_list, L, acceptDifferentReference = TRUE),
        tolerance = 1e-3)
})