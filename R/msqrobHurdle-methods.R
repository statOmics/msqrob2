
#' Function to fit msqrob hurdle models
#'
#' @noMd
#'
#' @description Fitting a hurdle msqrob model with an intensity component for
#'              assessing differential abundance and an count component that is
#'              modeling peptide counts using quasibinomial glm for differential
#'              detection of the number of features that were not missing and used
#'              for aggregation
#'
#' @param object \code{SummarizedExperiment} or \code{QFeatures} instance with an assay that is generated with the \code{aggreateFeatures} function from the \code{QFeatures} package
#'
#' @param formula Model formula. Both model components are built based on the
#'        covariates in the data object.
#'
#' @param modelColumnName \code{character} to indicate the variable name that is used
#'        to store the msqrob models in the rowData of the SummarizedExperiment
#'        instance or of the assay of the QFeatures instance. Default is "msqrobHurdle".
#'
#' @param overwrite \code{boolean(1)} to indicate if the column in the rowData has to
#'        be overwritten if the modelColumnName already exists. Default is FALSE.
#'
#' @param robust \code{boolean(1)} to indicate if robust regression is
#'     performed to account for outliers when fitting the intensity component of
#'     the hurdle model. Default is \code{TRUE}. If \code{FALSE} an OLS fit is performed.
#'
#' @param ridge \code{boolean(1)} to indicate if ridge regression is
#'        performed. Default is \code{FALSE}. If \code{TRUE} the fixed effects of the
#'        intensity component of the hurdle model are
#'        estimated via penalized regression and shrunken to zero.
#'
#' @param maxitRob \code{numeric(1)} indicating the maximum iterations in
#'        the IRWLS algorithm used in the M-estimation step of the robust
#'        regression for fitting the intensity component of the hurdle model.
#'
#' @param tol \code{numeric(1)} indicating the tolerance for declaring convergence
#'        of the M-estimation loop of the intensity component of the hurdle model.
#'
#' @param doQR \code{boolean(1)} to indicate if a QR decomposition is applied to the
#'     fixed-effect design matrix before the Scheipl ridge encoding. Default is \code{TRUE}.
#'     When \code{TRUE} the predictor columns are orthogonalised so that shrinkage is
#'     invariant to predictor ordering and collinearity. When \code{FALSE} the raw design
#'     matrix columns are used directly (standard L2 penalty on original parameters).
#'
#' @param lmerArgs a list (of correct class, resulting from ‘lmerControl()’
#'        containing control parameters, including the nonlinear optimizer to be used
#'        and parameters to be passed through to the nonlinear optimizer, see the
#'        ‘lmerControl’ documentation of the lme4 package for more details.
#'        Default is \code{list(control = lmerControl(calc.derivs = FALSE))}
#'
#' @param trend Character or FALSE. Controls the empirical Bayes variance trend
#'        for the intensity component. FALSE (default) or "none" uses a global
#'        prior variance. "mean" conditions on mean log-intensity. "count"
#'        conditions on log2 precursor count (from the .n column of rowData).
#'        "combined" uses a linear projection on both mean and log2 count.
#'
#' @param priorCount A 'numeric(1)', which is a prior count to be added to the observations to shrink
#'          the estimated odds ratios of the count component towards zero. Default is 0.1.
#'
#’ @param binomialBound logical, if ‘TRUE’ then the quasibinomial variance estimator will
#’                be never smaller than 1 (no underdispersion). Default is TRUE.
#’
#’ @examples
#'
#' # Load example data
#' # The data are a Feature object with containing
#' # a SummarizedExperiment named "peptide" with MaxQuant peptide intensities
#' # The data are a subset of spike-in the human-ecoli study
#' # The variable condition in the colData of the Feature object
#' # contains information on the spike in condition a-e (from low to high)
#' data(pe)
#'
#' # Aggregate peptide intensities to protein expression values
#' pe <- aggregateFeatures(pe, i = "peptide", fcol = "Proteins", name = "protein")
#'
#' # Fit Hurdle MSqrob model
#' # For summarized SummarizedExperiment
#' se <- pe[["protein"]]
#' se
#' colData(se) <- colData(pe)
#' se <- msqrobHurdle(se, formula = ~condition)
#' getCoef(rowData(se)$msqrobHurdleIntensity[[1]])
#' getCoef(rowData(se)$msqrobHurdleCount[[1]])
#'
#' # For features object
#' pe <- msqrobHurdle(pe, i = "protein", formula = ~condition)
#' getCoef(rowData(pe[["protein"]])$msqrobHurdleIntensity[[1]])
#' getCoef(rowData(pe[["protein"]])$msqrobHurdleCount[[1]])
#' @return SummarizedExperiment or QFeatures instance
#'
#' @rdname msqrobHurdle
#'
#' @aliases msqrobHurdle msqrobHurdle,SummarizedExperiment-method msqrobHurdle,QFeatures-method
#'
#' @author Lieven Clement
#'
#' @import SummarizedExperiment
#' @importFrom QFeatures aggcounts
#'
#' @export

setMethod(
    "msqrobHurdle", "SummarizedExperiment",
    function(object,
    formula,
    modelColumnName = "msqrobHurdle",
    overwrite = FALSE,
    robust = TRUE,
    ridge = FALSE,
    maxitRob = 1,
    tol = 1e-6,
    doQR = TRUE,
    trend = FALSE,
    lmerArgs = list(control = lmerControl(calc.derivs = FALSE)),
    priorCount = .1,
    binomialBound = TRUE) {
        if (ncol(colData(object)) == 0) stop("colData is empty")
        if ((sum(grep(pattern = modelColumnName, colnames(rowData(object)))) > 0) & !overwrite) {
            stop(
                "There are already columns with names starting with\'",
                modelColumnName,
                "\' in the rowData of the SummarizedExperiment object, set the argument overwrite=TRUE to replace the columns with the new results or use another name for the argument modelColumnName to store the results as novel columns in the rowData of SummarizedExperiment object"
            )
        }
        if (!(".n" %in% colnames(rowData(object)))) stop("The assay does not seem to be aggregated so the number of features used for aggregation are not available")
        if (!ridge) {
            rowData(object)[[paste0(modelColumnName, "Intensity")]] <- msqrobLm(
                y = assay(object),
                formula = formula,
                data = colData(object),
                robust = robust,
                maxitRob = maxitRob,
                trend = trend,
                counts = rowData(object)[[".n"]]
            )
        } else {
            rowData(object)[[paste0(modelColumnName, "Intensity")]] <- msqrobLmer(
                y = assay(object),
                formula = formula,
                data = colData(object),
                robust = robust,
                ridge = ridge,
                maxitRob = maxitRob,
                tol = tol,
                doQR = doQR,
                trend = trend,
                counts = rowData(object)[[".n"]],
                lmerArgs = lmerArgs
            )
        }
        rowData(object)[[paste0(modelColumnName, "Count")]] <- msqrobGlm(aggcounts(object),
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
#' @rdname msqrobHurdle

setMethod(
    "msqrobHurdle", "QFeatures",
    function(object,
    i,
    formula,
    modelColumnName = "msqrobHurdle",
    overwrite = FALSE,
    robust = TRUE,
    ridge = FALSE,
    maxitRob = 1,
    tol = 1e-6,
    doQR = TRUE,
    trend = FALSE,
    lmerArgs = list(control = lmerControl(calc.derivs = FALSE)),
    priorCount = .1,
    binomialBound = TRUE) {
        if (ncol(colData(object)) == 0) stop("colData is empty")
        if ((sum(grep(pattern = modelColumnName, colnames(rowData(object[[i]])))) > 0) & !overwrite) {
            stop(
                "There are already columns with names starting with\'",
                modelColumnName,
                "\' in the rowData of assay ",
                i,
                " of the QFeatures object, set the argument overwrite=TRUE to replace the columns with the new results or use another name for the argument modelColumnName to store the results as novel columns in the rowData"
            )
        }
        if (!(".n" %in% colnames(rowData(object[[i]])))) stop("The assay does not seem to be aggregated so the number of features used for aggregation is not available")
        if (!ridge) {
            rowData(object[[i]])[[paste0(modelColumnName, "Intensity")]] <- msqrobLm(
                y = assay(object[[i]]),
                formula = formula,
                data = colData(object),
                robust = robust,
                maxitRob = maxitRob,
                trend = trend,
                counts = rowData(object[[i]])[[".n"]]
            )
        } else {
            rowData(object[[i]])[[paste0(modelColumnName, "Intensity")]] <- msqrobLmer(
                y = assay(object[[i]]),
                formula = formula,
                data = colData(object),
                robust = robust,
                ridge = ridge,
                maxitRob = maxitRob,
                tol = tol,
                doQR = doQR,
                trend = trend,
                counts = rowData(object[[i]])[[".n"]],
                lmerArgs = lmerArgs
            )
        }
        rowData(object[[i]])[[paste0(modelColumnName, "Count")]] <- msqrobGlm(aggcounts(object[[i]]),
            rowData(object[[i]])[[".n"]],
            formula,
            colData(object),
            priorCount = priorCount,
            binomialBound = binomialBound
        )
        return(object)
    }
)
