#' @export
setGeneric("getModel", function(object) standardGeneric("getModel"))
#' @export
setGeneric("getFitMethod", function(object) standardGeneric("getFitMethod"))
#' @export
setGeneric("getCoef", function(object) standardGeneric("getCoef"))
#' @export
setGeneric("varContrast", function(object, ...) standardGeneric("varContrast"))
#' @export
setGeneric("getVar", function(object) standardGeneric("getVar"))
#' @export
setGeneric("getSigma", function(object) standardGeneric("getSigma"))
#' @export
setGeneric("getDF", function(object) standardGeneric("getDF"))
#' @export
setGeneric("getDfPosterior", function(object) standardGeneric("getDfPosterior"))
#' @export
setGeneric("getVarPosterior", function(object) standardGeneric("getVarPosterior"))
#' @export
setGeneric("getSigmaPosterior", function(object) standardGeneric("getSigmaPosterior"))
#' @export
setGeneric("getContrast", function(object, L) standardGeneric("getContrast"))
#' @export
setGeneric("getVcovUnscaled", function(object) standardGeneric("getVcovUnscaled"))
#' @export
setGeneric("msqrob", function(object, ...) standardGeneric("msqrob"))
#' @export
setGeneric("msqrobAggregate", function(object, ...) standardGeneric("msqrobAggregate"))
#' @export
setGeneric("hypothesisTest", function(object, ...) standardGeneric("hypothesisTest"))
#' @export
setGeneric("msqrobQB", function(object, ...) standardGeneric("msqrobQB"))
#' @export
setGeneric("msqrobHurdle", function(object, ...) standardGeneric("msqrobHurdle"))
#' @export
setGeneric("hypothesisTestHurdle", function(object, ...) standardGeneric("hypothesisTestHurdle"))
