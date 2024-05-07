#Create Data
#One protein that can be fit and one that returns a fitError
se_assay <- matrix(c(rnorm(n = 10,
                         mean = rep(c(5,10),
                                    each  =5),
                         sd = 0.1),
                   rep(NA,times=10)),
                   nrow =2,
                   byrow = TRUE,
                   dimnames = list(c("Protein1","Protein2"), LETTERS[1:10]))
se_cd <- DataFrame(treatment = rep(c("A","B"), each= 5),
                   row.names = LETTERS[1:10])
se <- SummarizedExperiment(assay = se_assay,
                           colData = se_cd)
qf <- QFeatures(List(se1 = se), colData = se_cd)
#Testing with SummarizedExperiment class
msqrob_output_se <- msqrob(object = se,
                           formula = ~ treatment,
                           robust = FALSE,
                           ridge = FALSE)

msqrob_output_qf <- msqrob(object = qf,
                           i = "se1",
                           formula = ~ treatment,
                           robust = FALSE,
                           ridge = FALSE)


test_that("Testing classes of msqrob output", {
  expect_s4_class(msqrob_output_se, "SummarizedExperiment")
  expect_s4_class(msqrob_output_qf, "QFeatures")
  expect_s4_class(msqrob_output_qf[["se1"]],"SummarizedExperiment")
})


#As the output of the SummarizedExperiment and QFeatures are made with
#The same functions i will not be testing them separately
#Todo: Testing all the functions separately

test_that("Testing msqrob models", {
  Protein1 <- rowData(msqrob_output_se)$msqrobModels$Protein1
  Protein2 <- rowData(msqrob_output_se)$msqrobModels$Protein2

  expect_s4_class(Protein1, "StatModel")
  expect_s4_class(Protein2, "StatModel")

  expect_equal(Protein1@type, expected = "lm")
  expect_equal(Protein2@type, expected = "fitError")
})



#Todo:
#Testing coldata
# - Missing data
# - Subsetted data
# - levels (ordering as well, important for coefficient names)
#Testing rowdata
#Testing parameters
# - correct naming (check ordering with levels of coldata!)
# - correct values
# ...
