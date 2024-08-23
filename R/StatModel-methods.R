#' Methods for StatModel class
#'
#' @description Methods for StatModel class
#'              \describe{
#'              \item{getContrast(object, L)}{to calculate contrasts of the model parameters}
#'              \item{varContrast(object, L)}{to calculate the variance-covariance matrix of the contrasts}
#'              }
#' @param object A list with elements of the class `StatModel` that are
#'        estimated using the \code{\link{msqrob}} function
#' @param L contrast `numeric` matrix specifying one or more contrasts of
#'        the linear model coefficients to be tested equal to zero.
#'        The rownames of the matrix should be equal to the names of
#'        parameters of the model.
#' @param acceptDifferentReference `boolean(1)` to indicate if the contrasts that involve
#'        parameters with modified reference levels are returned. Watch out putting this
#'        parameter to TRUE can change the interpretation of the logFC. Default is FALSE.
#'
#' @examples
#' data(pe)
#' # Aggregate peptide intensities in protein expression values
#' pe <- aggregateFeatures(pe, i = "peptide", fcol = "Proteins", name = "protein")
#'
#' # Fit msqrob model
#' pe <- msqrob(pe, i = "protein", formula = ~condition)
#'
#' # Define contrast
#' getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])
#' # Define contrast for log2 fold change between condition c and condition b:
#' L <- makeContrast("conditionc - conditionb=0", c("conditionb", "conditionc"))
#'
#' getContrast(rowData(pe[["protein"]])$msqrobModels[[1]], L)
#' varContrast(rowData(pe[["protein"]])$msqrobModels[[1]], L)
#' @rdname statModelMethods
#' @return A matrix with the calculated contrasts or variance-covariance matrix of contrasts
#' @aliases statModelMethods StatModel-method getContrast varContrast
#' @importFrom methods is

setMethod(
    "getContrast", "StatModel",
    function(object, L, acceptDifferentReference = FALSE) {
        if (!is(L, "matrix")) L <- as.matrix(L)
        if (referenceContrast(getReferencePresent(object), L, acceptDifferentReference)) {
            return(NA)
        }
        coefs <- getCoef(object)
        out <- matrix(rep(NA, ncol(L)))
        rownames(out) <- colnames(L)
        hlp <- try(t(L) %*% coefs[rownames(L)], silent = TRUE)
        if (!is(hlp, "try-error")) out[] <- hlp
        return(out)
    }
)

#' @rdname statModelMethods
#' @importFrom methods is

setMethod(
    "varContrast", "StatModel",
    function(object, L, acceptDifferentReference = FALSE) {
        if (!is(L, "matrix")) L <- as.matrix(L)
        if (referenceContrast(getReferencePresent(object), L, acceptDifferentReference)) {
            return(NA)
        }
        out <- matrix(NA, ncol(L), ncol(L))
        rownames(out) <- colnames(out) <- colnames(L)
        vcovTmp <- getVcovUnscaled(object) * object@varPosterior
        hlp <- try(t(L) %*% vcovTmp[rownames(L), rownames(L)] %*% L, silent = TRUE)
        if (!is(hlp, "try-error")) out[] <- hlp
        return(out)
    }
)
