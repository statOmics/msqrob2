setMethod("getModel",signature="msqrobModel",definition=function(object)
object@model)

setMethod("getFitMethod",signature="msqrobModel",definition=function(object)
object@modelType)

setMethod("getCoef",signature="msqrobModel",definition=function(object)
object@model$coefficients)

setMethod("getDF",signature="msqrobModel",definition=function(object)
if(object@modelType=="fitError") return(NA) else return(object@model$df.residual))
