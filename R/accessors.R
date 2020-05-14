#' Accessor functions for StatModel class
#'
#' @description Accessor functions for StatModel class
#'              \describe{
#'              \item{getModel(object)}{to get model}
#'              \item{getFitMethod(object)}{to get the parameter estimation method}
#'              \item{getCoef(object)}{to get the parameter estimates of the mean model}
#'              \item{getDF(object)}{to get the residual degrees of freedom of the model}
#'              \item{getVar(object)}{to get the residual variance of the model}
#'              \item{getSigma(object)}{to get the residual standard deviation of the model}
#'              \item{getDfPosterior(object)}{to get the degrees of freedom of
#'                the empirical Bayes variance estimator}
#'              \item{getVarPosterior(object)}{to get the empirical Bayes variance}
#'              \item{getSigmaPosterior(object)}{to get the empirical Bayes standard deviation}
#'              \item{getVcovUnscaled(object)}{to get the unscaled variance covariance matrix
#'                of the model parameters}
#'              }
#'
#' @rdname statModelAccessors
#' @aliases statModelAccessors getModel getFitMethod getCoef getDF getDfPosterior getVarPosterior getSigmaPosterior getVar getSigma getVcovUnscaled
#'
#' @param object `StatModel` object
setMethod("getModel",signature="StatModel",
          definition=function(object) object@params)

setMethod("getFitMethod",signature="StatModel",
          definition=function(object) object@type)

setMethod("getCoef",
          signature="StatModel",
          definition=function(object) object@params$coefficients)

setMethod("getDfPosterior",
          signature="StatModel",
          definition=function(object) object@dfPosterior)

setMethod("getVarPosterior",
          signature="StatModel",
          definition=function(object) object@varPosterior)

setMethod("getSigmaPosterior",
          signature="StatModel",
          definition=function(object) object@varPosterior^.5)

setMethod("getDF",
          signature="StatModel",
          definition=function(object) object@params$df.residual)

setMethod("getVar",
          signature="StatModel",
          definition=function(object)object@params$sigma^2)

setMethod("getSigma",
          signature="StatModel",
          definition=function(object)object@params$sigma)

setMethod("getVcovUnscaled",
          signature="StatModel",
          definition=function(object) object@params$vcovUnscaled)

#setMethod("getResults", "Features",
#          function(object, i, columnName){
#              if (!(i in names(object))) stop(paste0(i," is no assay of the Features object"))
#              if (is.null(columnName)) return(rowData(object[[i]]))
#              if (!(columnName %in% colnames(rowData(object[[i]])))) stop(paste0(columnName," is not a column of the rowData of the assay ",i," of the Features object"))
#              rowData(object[[i]])[[columnName]]
#              })
