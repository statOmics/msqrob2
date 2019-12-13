setMethod("getModel",signature="StatModel",definition=function(object)
object@params)

setMethod("getFitMethod",signature="StatModel",definition=function(object)
object@type)

setMethod("getCoef",signature="StatModel",definition=function(object)
object@params$coefficients)


setMethod("getDfPosterior",signature="StatModel",definition=function(object)
object@dfPosterior)

setMethod("getDF",signature="StatModel",definition=function(object)
{
if (object@type=="fitError") return(NA)
if (object@type=="rlm") return(sum(object@params$w)-object@params$rank)
else return(object@params$df.residual)
})
