
#' Function to fit msqrob models to peptide counts using glm
#'
#' @description Low-level function for parameter estimation with msqrob
#'              by modeling peptide counts using quasibinomial glm
#'
#' @param object `SummarizedExperiment` or `QFeatures` instance
#'
#' @param formula Model formula. The model is built based on the
#'        covariates in the data object.
#'
#' @param modelColumnName `character` to indicate the variable name that is used
#'        to store the msqrob models in the rowData of the SummarizedExperiment
#'        instance or of the assay of the QFeatures instance. Default is "msqrobModels".
#'
#' @param overwrite `boolean(1)` to indicate if the column in the rowData has to
#'        be overwritten if the modelColumnName already exists. Default is FALSE.
#'
#' @param priorCount A 'numeric(1)', which is a prior count to be added to the observations to shrink
#'          the estimated log-fold-changes towards zero. Default is 0.1.
#'
#' @param binomialBound logical, if ‘TRUE’ then the quasibinomial variance estimator will
#'                be never smaller than 1 (no underdispersion). Default is TRUE.
#'
#' @examples
#'
#' # Load example data
#' # The data are a Feature object with containing
#' # a SummarizedExperiment named "peptide" with MaxQuant peptide intensities
#' # The data are a subset of spike-in the human-ecoli study
#' # The variable condition in the colData of the Feature object
#' # contains information on the spike in condition a-e (from low to high)
#' data(pe)
#'
#' # Aggregate by counting how many peptide we observe for each protein
#' pe <- aggregateFeatures(pe, i = "peptide", fcol = "Proteins", name = "protein")
#'
#' # Fit MSqrob model to peptide counts using a quasi-binomial model
#' # For summarized SummarizedExperiment
#' se <- pe[["protein"]]
#' se
#' colData(se) <- colData(pe)
#' se <- msqrobQB(se, formula = ~condition)
#' getCoef(rowData(se)$msqrobQbModels[[1]])
#'
#' # For features object
#' pe <- msqrobQB(pe, i = "protein", formula = ~condition)
#' @return SummarizedExperiment or QFeatures instance
#'
#' @rdname msqrobQB
#'
#' @aliases msqrobQB msqrobQB,SummarizedExperiment-method msqrobQB,QFeatures-method
#'
#' @author Lieven Clement
#'
#' @import SummarizedExperiment
#' @importFrom QFeatures aggcounts
#'
#' @export

setMethod(
    "msqrobQB", "SummarizedExperiment",
    function(object,
    formula,
    modelColumnName = "msqrobQbModels",
    overwrite = FALSE,
    priorCount = .1,
    binomialBound = TRUE) {
        if (ncol(colData(object)) == 0) stop("colData is empty")
        if ((modelColumnName %in% colnames(rowData(object))) & !overwrite) {
            stop(
                "There is already a column named \'",
                modelColumnName,
                "\' in the rowData of the SummarizedExperiment object, set the argument overwrite=TRUE to replace the column with the new results or use another name for the argument modelColumnName to store the results as a novel column in the rowData of SummarizedExperiment object"
            )
        }
        if (!(".n" %in% colnames(rowData(object)))) stop("The assay does not seem to be aggregated so the number of features are not available")
        rowData(object)[[modelColumnName]] <- msqrobGlm(aggcounts(object),
            rowData(object)[[".n"]],
            formula,
            colData(object),
            priorCount = priorCount,
            binomialBound = binomialBound
        )
        return(object)
    }
)

#' @param i `character` or `integer` to specify the element of the `QFeatures` that
#'        contains the log expression intensities that will be modelled.
#' @export
#' @rdname msqrobQB

setMethod(
    "msqrobQB", "QFeatures",
    function(object,
    i,
    formula,
    modelColumnName = "msqrobQbModels",
    overwrite = FALSE,
    priorCount = .1,
    binomialBound = TRUE) {
        if (is.null(object[[i]])) stop("QFeatures object does not contain an assay with the name ", i)
        if ((modelColumnName %in% colnames(rowData(object[[i]]))) & !overwrite) {
            stop(
                "There is already a column named \'",
                modelColumnName,
                "\' in the rowData of assay \'",
                i,
                "'of object, set the argument overwrite=TRUE to replace the column with the new results or use another name for the argument modelColumnName to store the results as a novel column in the rowData of SummarizedExperiment object"
            )
        }
        if (!(".n" %in% colnames(rowData(object[[i]])))) stop("The assay does not seem to be aggregated so the number of features are not available")
        rowData(object[[i]])[[modelColumnName]] <- msqrobGlm(aggcounts(object[[i]]),
            rowData(object[[i]])[[".n"]],
            formula,
            colData(object),
            priorCount = priorCount,
            binomialBound = binomialBound
        )
        return(object)
    }
)
