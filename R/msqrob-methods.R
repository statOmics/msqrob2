#' Methods to fit msqrob models with ridge regression and/or random effects using lme4
#'
#' @description Parameter estimation of msqrob models for `QFeatures`
#'              and `SummarizedExperiment` instance.
#'
#' @aliases msqrob msqrob,SummarizedExperiment-method msqrob,QFeatures-method
#'
#' @author Lieven Clement
#'
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
#' # Aggregate peptide intensities in protein expression values
#' pe <- aggregateFeatures(pe, i = "peptide", fcol = "Proteins", name = "protein")
#'
#' # Fit MSqrob model using robust linear regression upon summarization of
#' # peptide intensities into protein expression values.
#' # For summarized SummarizedExperiment
#' se <- pe[["protein"]]
#' se
#' colData(se) <- colData(pe)
#' se <- msqrob(se, formula = ~condition, modelColumnName = "rlm")
#' getCoef(rowData(se)$rlm[[1]])
#'
#' # For features object
#' pe <- msqrob(pe, i = "protein", formula = ~condition, modelColumnName = "rlm")
#' # with ridge regression (slower)
#' pe <- msqrob(pe, i = "protein", formula = ~condition, ridge = TRUE, modelColumnName = "ridge")
#'
#' # compare for human protein (no DE)==> large shrinkage to zero
#' cbind(getCoef(rowData(pe[["protein"]])$rlm[[1]]), getCoef(rowData(pe[["protein"]])$ridge[[1]]))
#'
#' # compare for ecoli protein (DE)==> almost no shrinkage to zero
#' cbind(
#'     getCoef(rowData(pe[["protein"]])$rlm[["P00956"]]),
#'     getCoef(rowData(pe[["protein"]])$ridge[["P00956"]])
#' )
#' @param object `SummarizedExperiment` or `QFeatures` instance
#'
#' @param formula Model formula. The model is built based on the
#'     covariates in the data object.
#'
#' @param modelColumnName `character` to indicate the variable name that is used
#'        to store the msqrob models in the rowData of the SummarizedExperiment
#'        instance or of the assay of the QFeatures instance. Default is "msqrobModels".
#'
#' @param overwrite `boolean(1)` to indicate if the column in the rowData has to
#'        be overwritten if the modelColumnName already exists. Default is FALSE.
#'
#' @param robust `boolean(1)` to indicate if robust regression is
#'     performed to account for outliers. Default is `TRUE`. If
#'     `FALSE` an OLS fit is performed.
#'
#' @param ridge `boolean(1)` to indicate if ridge regression is
#'        performed. Default is `FALSE`. If `TRUE` the fixed effects are
#'        estimated via penalized regression and shrunken to zero.
#'
#' @param maxitRob `numeric(1)` indicating the maximum iterations in
#'        the IRWLS algorithm used in the M-estimation step of the robust
#'        regression.
#'
#' @param tol `numeric(1)` indicating the tolerance for declaring convergence
#'        of the M-estimation loop.
#'
#' @param doQR `boolean(1)` to indicate if QR decomposition is used when adopting
#'     ridge regression. Default is `TRUE`. If `FALSE` the predictors of the fixed
#'     effects are not transformed, and the degree of shrinkage can depend on the encoding.
#'
#' @param lmerArgs a list (of correct class, resulting from ‘lmerControl()’
#'        containing control parameters, including the nonlinear optimizer to be used
#'        and parameters to be passed through to the nonlinear optimizer, see the
#'        ‘lmerControl’ documentation of the lme4 package for more details.
#'        Default is `list(control = lmerControl(calc.derivs = FALSE))`
#' @rdname msqrob
#'
#' @import SummarizedExperiment
#' @export
setMethod(
    "msqrob", "SummarizedExperiment",
    function(object,
    formula,
    modelColumnName = "msqrobModels",
    overwrite = FALSE,
    robust = TRUE,
    ridge = FALSE,
    maxitRob = 1,
    tol = 1e-6,
    doQR = TRUE,
    lmerArgs = list(control = lmerControl(calc.derivs = FALSE))) {
        if (ncol(colData(object)) == 0) stop("colData is empty")
        if ((modelColumnName %in% colnames(rowData(object))) & !overwrite) {
            stop(
                "There is already a column named \'",
                modelColumnName,
                "\' in the rowData of the SummarizedExperiment object, set the argument overwrite=TRUE to replace the column with the new results or use another name for the argument modelColumnName to store the results as a novel column in the rowData of SummarizedExperiment object"
            )
        }
        if (!ridge) {
            rowData(object)[[modelColumnName]] <- msqrobLm(
                y = assay(object),
                formula = formula,
                data = colData(object),
                robust = robust,
                maxitRob = maxitRob
            )
        } else {
            rowData(object)[[modelColumnName]] <- msqrobLmer(
                y = assay(object),
                formula = formula,
                data = colData(object),
                robust = robust,
                maxitRob = maxitRob,
                tol = tol,
                doQR = doQR,
                lmerArgs = lmerArgs
            )
        }
        return(object)
    }
)

#' @param i `character` or `integer` to specify the element of the `QFeatures` that
#'        contains the log expression intensities that will be modelled.
#'
#' @return A SummarizedExperiment or a `QFeatures` instance with the models.
#' @export
#' @rdname msqrob
setMethod(
    "msqrob", "QFeatures",
    function(object,
    i,
    formula,
    modelColumnName = "msqrobModels",
    overwrite = FALSE,
    robust = TRUE,
    ridge = FALSE,
    maxitRob = 1,
    tol = 1e-6,
    doQR = TRUE,
    lmerArgs = list(control = lmerControl(calc.derivs = FALSE))) {
        if (is.null(object[[i]])) stop("QFeatures object does not contain an assay with the name ", i)
        if ((modelColumnName %in% colnames(rowData(object[[i]]))) & !overwrite) {
            stop(
                "There is already a column named \'",
                modelColumnName,
                "\' in the rowData of the assay",
                i,
                "of the QFeatures object, set the argument overwrite=TRUE to replace the column with the new results or use another name for the argument modelColumnName to store the results as a novel column in the rowData of assay of the QFeatures object"
            )
        }
        if (!ridge) {
            rowData(object[[i]])[[modelColumnName]] <- msqrobLm(
                y = assay(object[[i]]),
                formula = formula,
                data = colData(object),
                robust = robust,
                maxitRob = maxitRob
            )
        } else {
            rowData(object[[i]])[[modelColumnName]] <- msqrobLmer(
                y = assay(object[[i]]),
                formula = formula,
                data = colData(object),
                robust = robust,
                maxitRob = maxitRob,
                tol = tol,
                doQR = doQR,
                lmerArgs = lmerArgs
            )
        }
        return(object)
    }
)
