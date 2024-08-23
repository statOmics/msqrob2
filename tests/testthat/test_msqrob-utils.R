.create_minimal_data <- function() {
    y_no_ref <- c(rep(NA,5), rep(1,10))
    names(y_no_ref) <- 1:15
    y_ref <- c(rep(1,5), rep(1,5), rep(1,5))
    names(y_ref) <- 1:15
    data_no_ref <- data.frame(condition = as.factor(rep(letters[1:3], c(5,5,5))),
        numerical = c(1:5, 1:5, 1:5),
        row.names = 1:15)
    form <- formula(~ 1 + condition + numerical, data = data)
    form_num <- formula(~ 1 + numerical, data = data)
    return(list(data = data, form = form,
        form_num = form_num,
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

    referenceCond <- "a"
    names(referenceCond) <- "condition"

    reference_present_no_ref <- c(FALSE, FALSE)
    names(reference_present_no_ref) <- c("conditionb", "conditionc")
    expect_identical(reference_present_no_ref, msqrob2:::checkReference(y_no_ref, data, referenceCond))

    reference_present_ref <- c(TRUE, TRUE)
    names(reference_present_ref) <- c("conditionb", "conditionc")
    expect_identical(reference_present_ref, msqrob2:::checkReference(y_ref, data, referenceCond))

    expect_identical(NULL, msqrob2:::checkReference(data, form, NULL))
})

test_that("referenceContrast", {
    L <- matrix(c(1, 0), nrow = 2, byrow = TRUE)
    colnames(L) <- "conditionb=0"
    rownames(L) <- c("conditionb", "conditionc")

    reference_present_no_ref <- c(FALSE, FALSE)
    names(reference_present_no_ref) <- c("conditionb", "conditionc")

    expect_identical(TRUE, msqrob2:::referenceContrast(reference_present_no_ref, L, FALSE))

    reference_present_ref <- c(TRUE, TRUE)
    names(reference_present_ref) <- c("conditionb", "conditionc")

    expect_identical(FALSE, msqrob2:::referenceContrast(reference_present_ref, L, FALSE))

    expect_identical(FALSE, msqrob2:::referenceContrast(reference_present_no_ref, L, TRUE))
    expect_identical(FALSE, msqrob2:::referenceContrast(reference_present_ref, L, TRUE))
})