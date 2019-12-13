setMethod("getModel",signature="msqrobModel",definition=function(object)
object@model)

setMethod("getFitMethod",signature="msqrobModel",definition=function(object)
object@modelType)

setMethod("getCoef",signature="msqrobModel",definition=function(object)
object@model$coefficients)


setMethod("getDfPosterior",signature="msqrobModel",definition=function(object)
object@dfPosterior)

setMethod("getDF",signature="msqrobModel",definition=function(object)
{
if (object@modelType=="fitError") return(NA)
if (object@modelType=="rlm") return(sum(object@model$w)-object@model$rank)
else return(object@model$df.residual)
})
