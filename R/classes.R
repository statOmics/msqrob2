#' The msqrobModel class
#'
#' This class represents a model fitted to a feature in a proteomics experiment using msqrob
#'
#' @section Slots:
#' \describe{
#'  \item{\code{modelType}:}{msqrob provides four options to estimate a model. ols: estimation via ordinary least squares, robust: parameter estimation via robust regression, ridge: penalised estimation using a ridge penalty, robustRidge: penalised M-estimation using a ridge penalty.}
#' \item{\code{model}:}{Model object. A list containing the estimated model parameters required for downstream inference.}
#' }
#' @name msqrobModel-class
#' @rdname msqrobModel-class
#' @exportClass msqrobModel
#'

setClass("msqrobModel",representation=representation(modelType="character",model="list"),prototype = prototype(modelType="",model=list()))
