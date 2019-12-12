setMethod("getModel",signature="msqrobModel",definition=function(object)
object@model)

setMethod("getFitMethod",signature="msqrobModel",definition=function(object)
object@modelType)

setMethod("getCoef",signature="msqrobModel",definition=function(object)
object@model$coefficients)
