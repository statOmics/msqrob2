.create_minimal_data <- function() {
    y_no_ref <- c(rep(NA,5), rep(1,10))
    names(y_no_ref) <- 1:15
    y_ref <- c(rep(1,5), rep(1,5), rep(1,5))
    names(y_ref) <- 1:15
    data <- data.frame(condition = as.factor(rep(letters[1:3], c(5,5,5))),
        condition2 = as.factor(rep(letters[1:5], 3)),
        numerical = c(1:5, 1:5, 1:5),
        row.names = 1:15)
    form <- formula(~ 1 + condition + numerical, data = data)
    form_num <- formula(~ 1 + numerical, data = data)
    form_no_intercept <- formula(~ -1 + condition + numerical, data = data)
    form_interaction <- formula(~ 1 + condition*condition2, data = data)
    return(list(data = data, form = form,
        form_num = form_num,
        form_no_intercept = form_no_intercept,
        form_interaction = form_interaction,
        y_no_ref = y_no_ref,
        y_ref = y_ref))
}

test_that("getReferenceLevels", {
    data <- .create_minimal_data()$data
    form <- .create_minimal_data()$form
    form_num <- .create_minimal_data()$form_num

    expect_identical(NULL, msqrob2:::getReferenceLevels(data, form_num))

    referenceCond <- "a"
    names(referenceCond) <- "condition"
    expect_identical(referenceCond, msqrob2:::getReferenceLevels(data, form))
})

test_that("checkReference", {
    data <- .create_minimal_data()$data
    y_no_ref <- .create_minimal_data()$y_no_ref
    y_ref <- .create_minimal_data()$y_ref
    formula <- .create_minimal_data()$form
    formula_no_intercept <- .create_minimal_data()$form_no_intercept
    formula_interaction <- .create_minimal_data()$form_interaction

    paramNames <- colnames(model.matrix(formula, data = data))
    paramNames_no_intercept <- colnames(model.matrix(formula_no_intercept, data = data))
    paramNames_interaction <- colnames(model.matrix(formula_interaction, data = data))

    reference_present_no_ref <- c(FALSE, FALSE)
    names(reference_present_no_ref) <- c("conditionb", "conditionc")
    expect_identical(reference_present_no_ref, msqrob2:::checkReference(y_no_ref, data, paramNames, formula))

    reference_present_no_int <- c(TRUE, TRUE, TRUE)
    names(reference_present_no_int) <- c("conditiona", "conditionb", "conditionc")
    expect_identical(reference_present_no_int, msqrob2:::checkReference(y_no_ref, data, paramNames_no_intercept, formula_no_intercept))

    reference_present_ref <- c(TRUE, TRUE)
    names(reference_present_ref) <- c("conditionb", "conditionc")
    expect_identical(reference_present_ref, msqrob2:::checkReference(y_ref, data, paramNames, formula))

    reference_no_var <- logical(0)
    expect_identical(reference_no_var, msqrob2:::checkReference(y_ref, data, c("(Intercept)"), as.formula(~1)))

    reference_interaction <- rep(FALSE, length(paramNames_interaction[-1]))
    names(reference_interaction) <- paramNames_interaction[-1]
    expect_identical(reference_interaction, msqrob2:::checkReference(y_no_ref, data, paramNames_interaction, formula_interaction))
})

test_that("referenceContrast", {
    L <- matrix(c(1, 0), nrow = 2, byrow = TRUE)
    colnames(L) <- "conditionb=0"
    rownames(L) <- c("conditionb", "conditionc")

    reference_present_no_ref <- c(FALSE, FALSE)
    names(reference_present_no_ref) <- c("conditionb", "conditionc")

    expect_identical(FALSE, msqrob2:::referenceContrast(reference_present_no_ref, L, FALSE))

    reference_present_ref <- c(TRUE, TRUE)
    names(reference_present_ref) <- c("conditionb", "conditionc")

    expect_identical(TRUE, msqrob2:::referenceContrast(reference_present_ref, L, FALSE))

    expect_identical(TRUE, msqrob2:::referenceContrast(reference_present_no_ref, L, TRUE))
    expect_identical(TRUE, msqrob2:::referenceContrast(reference_present_ref, L, TRUE))
})