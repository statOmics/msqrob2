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
#' @aliases statModelAccessors getCoef getDF getDfPosterior getFitMethod getModel getSigma getSigmaPosterior getVar getVarPosterior getVcovUnscaled
#'
#' @param object `StatModel` object

setMethod("getModel",
    signature = "StatModel",
    definition = function(object) object@params
)


#' @rdname statModelAccessors
setMethod("getFitMethod",
    signature = "StatModel",
    definition = function(object) object@type
)

#' @rdname statModelAccessors
setMethod("getCoef",
    signature = "StatModel",
    definition = function(object) object@params$coefficients
)

#' @rdname statModelAccessors
setMethod("getDfPosterior",
    signature = "StatModel",
    definition = function(object) object@dfPosterior
)

#' @rdname statModelAccessors
setMethod("getVarPosterior",
    signature = "StatModel",
    definition = function(object) object@varPosterior
)

#' @rdname statModelAccessors
setMethod("getSigmaPosterior",
    signature = "StatModel",
    definition = function(object) object@varPosterior^.5
)

#' @rdname statModelAccessors
setMethod("getDF",
    signature = "StatModel",
    definition = function(object) object@params$df.residual
)

#' @rdname statModelAccessors
setMethod("getVar",
    signature = "StatModel",
    definition = function(object) object@params$sigma^2
)

#' @rdname statModelAccessors
setMethod("getSigma",
    signature = "StatModel",
    definition = function(object) object@params$sigma
)

#' @rdname statModelAccessors
setMethod("getVcovUnscaled",
    signature = "StatModel",
    definition = function(object) object@params$vcovUnscaled
)
