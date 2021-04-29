#' Parameter estimates, standard errors and statistical inference on differential
#' expression analysis
#'
#' @description Summary table of the estimates for differential expression of features
#'
#' @rdname hypothesisTest
#'
#' @aliases hypothesisTest hypothesisTest,SummarizedExperiment-method hypothesisTest,QFeatures-method hypothesisTestHurdle hypothesisTestHurdle,SummarizedExperiment-method hypothesisTestHurdle,QFeatures-method
#'
#' @author Lieven Clement
#'
#' @examples
#'
#' # Load example data
#' # The data are a Feature object containing
#' # a SummarizedExperiment named "peptide" with MaxQuant peptide intensities
#' # The data are a subset of spike-in the human-ecoli study
#' # The variable condition in the colData of the Feature object
#' # contains information on the spike in condition a-e (from low to high)
#' data(pe)
#'
#' # Aggregate peptide intensities in protein expression values
#' pe <- aggregateFeatures(pe, i = "peptide", fcol = "Proteins", name = "protein")
#'
#' # Fit msqrob model
#' pe <- msqrob(pe, i = "protein", formula = ~condition)
#'
#' # Define contrast
#' getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])
#' # Assess log2 fold change between condition c and condition b
#' L <- makeContrast(
#'     "conditionc - conditionb=0",
#'     c("conditionb", "conditionc")
#' )
#'
#' # example SummarizedExperiment instance
#' se <- pe[["protein"]]
#' se <- hypothesisTest(se, L)
#' head(rowData(se)$"conditionc - conditionb", 10)
#' # Volcano plot
#' plot(-log10(pval) ~ logFC,
#'     rowData(se)$"conditionc - conditionb",
#'     col = (adjPval < 0.05) + 1
#' )
#'
#' # Example for QFeatures instance
#' # Assess log2 fold change between condition b and condition a (reference class),
#' # condition c and condition a, and, condition c and condition b.
#' L <- makeContrast(
#'     c(
#'         "conditionb=0",
#'         "conditionc=0",
#'         "conditionc - conditionb=0"
#'     ),
#'     c("conditionb", "conditionc")
#' )
#' pe <- hypothesisTest(pe, i = "protein", L)
#' head(rowData(pe[["protein"]])$"conditionb", 10)
#' # Volcano plots
#' par(mfrow = c(1, 3))
#' plot(-log10(pval) ~ logFC,
#'     rowData(pe[["protein"]])$"conditionb",
#'     col = (adjPval < 0.05) + 1,
#'     main = "log2 FC b-a"
#' )
#' plot(-log10(pval) ~ logFC,
#'     rowData(pe[["protein"]])$"conditionc",
#'     col = (adjPval < 0.05) + 1,
#'     main = "log2 FC c-a"
#' )
#' plot(-log10(pval) ~ logFC,
#'     rowData(pe[["protein"]])$"conditionc - conditionb",
#'     col = (adjPval < 0.05) + 1,
#'     main = "log2 FC c-b"
#' )
#'
#' # Hurdle method
#' pe <- msqrobHurdle(pe, i = "protein", formula = ~condition)
#' pe <- hypothesisTestHurdle(pe, i = "protein", L)
#' head(rowData(pe[["protein"]])$"hurdle_conditionb", 10)
#' @param object `SummarizedExperiment` or `QFeatures` instance
#' @param contrast `numeric` matrix specifying one or more contrasts of
#'        the linear model coefficients to be tested equal to zero. If multiple
#'        contrasts are given (multiple columns) then results will be returned for
#'        each contrast. The rownames of the matrix should be equal to the names
#'        of parameters of the model that are involved in the contrast.
#'        The column names of the matrix will be used to construct names to store
#'        the results in the rowData of the SummarizedExperiment or of the assay of
#'        the QFeatures object. The contrast matrix can be made using the `makeContrast`
#'        function.
#' @param adjust.method `character` specifying the method to adjust
#'        the p-values for multiple testing.
#'        Options, in increasing conservatism, include ‘"none"’,
#'        ‘"BH"’, ‘"BY"’ and ‘"holm"’.  See ‘p.adjust’ for the complete
#'        list of options. Default is "BH" the Benjamini-Hochberg method
#'        to controle the False Discovery Rate (FDR).
#' @param modelColumn `character` to indicate the variable name that was used
#'        to store the msqrob models in the rowData of the SummarizedExperiment
#'        instance or of the assay of the QFeatures instance. Default is "msqrobModels"
#'        when the `hypothesisTest` function is used and "msqrobHurdle" for `hypothesisTestHurdle`.
#' @param resultsColumnNamePrefix `character` to indicate the the prefix for the
#'        variable name that will be used to store test results in the rowData of
#'        the SummarizedExperiment instance or of the assay of the QFeatures instance.
#'        Default is "" so that the variable name with the results will be
#'        the column name of the column in the contrast matrix L. If L is a matrix
#'        with multiple columns, multiple results columns will be made, one for each
#'        contrast. If L is a matrix with a single column which has no column names and if resultsColumnNamePrefix="" the
#'        results will be stored in the column with name msqrobResults. For hypothesisTestHurdle
#'        the default prefix is "hurdle_". If L is a matrix with one column and has no column names and if resultsColumnNamePrefix="hurdle_" the
#'        results will be stored in the column with name hurdleResults.
#' @param overwrite `boolean(1)` to indicate if the column in the rowData has to
#'        be overwritten if the modelColumnName already exists. Default is FALSE.
#'
#' @import SummarizedExperiment
#' @importFrom stats pchisq p.adjust
#' @export

setMethod(
    "hypothesisTest", "SummarizedExperiment",
    function(object,
    contrast,
    adjust.method = "BH",
    modelColumn = "msqrobModels",
    resultsColumnNamePrefix = "",
    overwrite = FALSE) {
        if (!(modelColumn %in% colnames(rowData(object)))) stop("There is no column named \'", modelColumn, "\' with stored models of an msqrob fit in the rowData of the SummarizedExperiment object")
        if (is.null(colnames(contrast)) & resultsColumnNamePrefix == "") resultsColumnNamePrefix <- "msqrobResults"
        if (is.null(colnames(contrast)) & ncol(contrast) > 1) colnames(contrast) <- seq_len(ncol(contrast))
        if ((sum(paste0(resultsColumnNamePrefix, colnames(contrast)) %in% colnames(rowData(object))) > 0) & !overwrite) stop("There is/are already column(s) with names starting with\'", resultsColumnNamePrefix, "\' in the rowData of the SummarizedExperiment object, set the argument overwrite=TRUE to replace the column(s) with the new results or use another name for the argument resultsColumnNamePrefix")
        for (j in seq_len(ncol(contrast)))
        {
            contrHlp <- contrast[, j]
            names(contrHlp) <- rownames(contrast)
            rowData(object)[[paste0(resultsColumnNamePrefix, colnames(contrast)[j])]] <- topFeatures(rowData(object)[, modelColumn], contrast = contrHlp, adjust.method = adjust.method, sort = FALSE, alpha = 1)
        }
        return(object)
    }
)

#' @import SummarizedExperiment
#' @importFrom stats pchisq p.adjust
#' @export
#' @rdname hypothesisTest

setMethod(
    "hypothesisTestHurdle", "SummarizedExperiment",
    function(object,
    contrast,
    adjust.method = "BH",
    modelColumn = "msqrobHurdle",
    resultsColumnNamePrefix = "hurdle_",
    overwrite = FALSE) {
        if (sum(paste0(modelColumn, c("Intensity", "Count")) %in% colnames(rowData(object))) != 2) stop("There are no columns for the models of the hurdle components in the rowData of the SummarizedExperiment")
        if (is.null(colnames(contrast)) & resultsColumnNamePrefix == "hurdle_") resultsColumnNamePrefix <- "hurdleResults"
        if (is.null(colnames(contrast)) & ncol(contrast) > 1) colnames(contrast) <- seq_len(ncol(contrast))
        if ((sum(paste0(resultsColumnNamePrefix, colnames(contrast)) %in% colnames(rowData(object))) > 0) & !overwrite) stop("There is/are already column(s) with names starting with\'", resultsColumnNamePrefix, "\' in the rowData of the SummarizedExperiment object, set the argument overwrite=TRUE to replace the column(s) with the new results or use another name for the argument resultsColumnNamePrefix")
        for (j in seq_len(ncol(contrast)))
        {
            contrHlp <- contrast[, j]
            names(contrHlp) <- rownames(contrast)
            intensityComponent <- topFeatures(rowData(object)[, paste0(modelColumn, "Intensity")], contrast = contrHlp, adjust.method = adjust.method, sort = FALSE, alpha = 1)
            countComponent <- topFeatures(rowData(object)[, paste0(modelColumn, "Count")], contrast = contrHlp, adjust.method = adjust.method, sort = FALSE, alpha = 1)

            sam <- cbind(
                intensityComponent[, seq_len(5)],
                countComponent[, seq_len(5)]
            )

            colnames(sam)[seq(2, 5)] <- paste0("logFC", colnames(sam)[seq(2, 5)])
            colnames(sam)[6] <- "logOR"
            colnames(sam)[seq(7, 10)] <- paste0("logOR", colnames(sam)[seq(7, 10)])

            sam$fisher <- -2 * (log(sam[, 5]) + log(sam[, 10]))
            sam$fisherDf <- 4
            sam$fisherDf[is.na(sam$fisher)] <- 2
            id1 <- is.na(sam$fisher) & !is.na(sam[, 5])
            id2 <- is.na(sam$fisher) & !is.na(sam[, 10])
            sam$fisher[id1] <- -2 * log(sam[id1, 5])
            sam$fisher[id2] <- -2 * log(sam[id2, 10])

            sam$fisherPval <- pchisq(sam$fisher, sam$fisherDf, lower.tail = FALSE)
            sam$fisherAdjPval <- p.adjust(sam$fisherPval, adjust.method)

            rowData(object)[[paste0(resultsColumnNamePrefix, colnames(contrast)[j])]] <- sam
        }
        return(object)
    }
)


#' @param i `character` or `integer` to specify the element of the `QFeatures` that
#'        contains the log expression intensities that will be modelled.
#'
#' @return A SummarizedExperiment or a `QFeatures` instance augmented with the test
#'         results.
#'
#' @export
#' @rdname hypothesisTest

setMethod(
    "hypothesisTest", "QFeatures",
    function(object,
    i,
    contrast,
    adjust.method = "BH",
    modelColumn = "msqrobModels",
    resultsColumnNamePrefix = "",
    overwrite = FALSE) {
        if (is.null(object[[i]])) stop("QFeatures object does not contain an assay with the name ", i)
        if (!(modelColumn %in% colnames(rowData(object[[i]])))) stop("There is no column named \'", modelColumn, "\' with stored models of an msqrob fit in the rowData of assay ", i, "of the QFeatures object.")
        if (is.null(colnames(contrast)) & resultsColumnNamePrefix == "") resultsColumnNamePrefix <- "msqrobResults"
        if (is.null(colnames(contrast)) & ncol(contrast) > 1) colnames(contrast) <- seq_len(ncol(contrast))
        if ((sum(paste0(resultsColumnNamePrefix, colnames(contrast)) %in% colnames(rowData(object[[i]]))) > 0) & !overwrite) stop("There is/are already column(s) with names starting with", resultsColumnNamePrefix, "\' in the rowData of assay ", i, " of the QFeatures object, set the argument overwrite=TRUE to replace the column(s) with the new results or use another name for the argument resultsColumnNamePrefix")
        for (j in seq_len(ncol(contrast)))
        {
            contrHlp <- contrast[, j]
            names(contrHlp) <- rownames(contrast)
            rowData(object[[i]])[[paste0(resultsColumnNamePrefix, colnames(contrast)[j])]] <- topFeatures(rowData(object[[i]])[, modelColumn], contrast = contrHlp, adjust.method = adjust.method, sort = FALSE, alpha = 1)
        }
        return(object)
    }
)

#' @export
#' @rdname hypothesisTest

setMethod(
    "hypothesisTestHurdle", "QFeatures",
    function(object,
    i,
    contrast,
    adjust.method = "BH",
    modelColumn = "msqrobHurdle",
    resultsColumnNamePrefix = "hurdle_",
    overwrite = FALSE) {
        if (is.null(object[[i]])) stop("QFeatures object does not contain an assay with the name ", i)
        if (sum(paste0(modelColumn, c("Intensity", "Count")) %in% colnames(rowData(object[[i]]))) != 2) stop("There are no columns for the models of the hurdle components in the rowData of assay ", i, "of the QFeatures object.")
        if (is.null(colnames(contrast)) & resultsColumnNamePrefix == "hurdle_") resultsColumnNamePrefix <- "hurdleResults"
        if (is.null(colnames(contrast)) & ncol(contrast) > 1) colnames(contrast) <- seq_len(ncol(contrast))
        if ((sum(paste0(resultsColumnNamePrefix, colnames(contrast)) %in% colnames(rowData(object[[i]]))) > 0) & !overwrite) stop("There is/are already column(s) with names starting with ", resultsColumnNamePrefix, "\' in the rowData of assay ", i, " of the QFeatures object, set the argument overwrite=TRUE to replace the column(s) with the new results or use another name for the argument resultsColumnNamePrefix")
        for (j in seq_len(ncol(contrast)))
        {
            contrHlp <- contrast[, j]
            names(contrHlp) <- rownames(contrast)
            intensityComponent <- topFeatures(rowData(object[[i]])[, paste0(modelColumn, "Intensity")], contrast = contrHlp, adjust.method = adjust.method, sort = FALSE, alpha = 1)
            countComponent <- topFeatures(rowData(object[[i]])[, paste0(modelColumn, "Count")], contrast = contrHlp, adjust.method = adjust.method, sort = FALSE, alpha = 1)

            sam <- cbind(
                intensityComponent[, seq_len(5)],
                countComponent[, seq_len(5)]
            )

            colnames(sam)[seq(2, 5)] <- paste0("logFC", colnames(sam)[seq(2, 5)])
            colnames(sam)[6] <- "logOR"
            colnames(sam)[seq(7, 10)] <- paste0("logOR", colnames(sam)[seq(7, 10)])

            sam$fisher <- -2 * (log(sam[, 5]) + log(sam[, 10]))
            sam$fisherDf <- 4
            sam$fisherDf[is.na(sam$fisher)] <- 2
            id1 <- is.na(sam$fisher) & !is.na(sam[, 5])
            id2 <- is.na(sam$fisher) & !is.na(sam[, 10])
            sam$fisher[id1] <- -2 * log(sam[id1, 5])
            sam$fisher[id2] <- -2 * log(sam[id2, 10])

            sam$fisherPval <- pchisq(sam$fisher, sam$fisherDf, lower.tail = FALSE)
            sam$fisherAdjPval <- p.adjust(sam$fisherPval, adjust.method)

            rowData(object[[i]])[[paste0(resultsColumnNamePrefix, colnames(contrast)[j])]] <- sam
        }
        return(object)
    }
)
