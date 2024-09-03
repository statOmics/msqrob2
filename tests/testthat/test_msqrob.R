.create_minimal_data <- function() {
    y_no_ref <- c(rep(NA,5), runif(10))
    y_ref <- c(runif(5), runif(5), runif(5))
    Y <- matrix(c(y_no_ref, y_ref), nrow = 2, ncol = 15, byrow = TRUE)
    Y_defficient <- matrix(c(y_no_ref, y_no_ref), nrow = 2, ncol = 15, byrow = TRUE)
    rownames(Y) <- c("feat1", "feat2")
    colnames(Y) <- paste0("S", 1:15)
    data <- data.frame(condition = as.factor(rep(letters[1:3], c(5,5,5))),
        numerical = c(1:5, 1:5, 1:5),
        row.names = paste0("S", 1:15))
    form_cond <- formula(~ 1 + condition, data = data)
    return(list(data = data, form_cond = form_cond, Y = Y, Y_defficient = Y_defficient))
}

test_that("msqrobLm", {
    set.seed(123)
    Y_defficient <- .create_minimal_data()$Y_defficient
    data <- .create_minimal_data()$data
    form_cond <- .create_minimal_data()$form_cond
    expect_error(msqrobLm(Y_defficient, form_cond, data, robust = FALSE),
        "The provided model must be full rank")

    Y <- .create_minimal_data()$Y
    msqrobLm_object <- msqrobLm(Y, form_cond, data, robust = FALSE)

    reference_present_no_ref <- c(FALSE, FALSE)
    names(reference_present_no_ref) <- c("conditionb", "conditionc")
    expect_identical(reference_present_no_ref, msqrobLm_object$feat1@params$referencePresent)

    reference_present_ref <- c(TRUE, TRUE)
    names(reference_present_ref) <- c("conditionb", "conditionc")
    expect_identical(reference_present_ref, msqrobLm_object$feat2@params$referencePresent)
    })
