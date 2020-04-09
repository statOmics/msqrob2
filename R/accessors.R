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
return(object@params$df.residual)
})

setMethod("getVar",signature="StatModel",definition=function(object)
{
return(object@params$sigma^2)
})

setMethod("getVcovUnscaled",signature="StatModel",definition=function(object)
{
return(object@params$vcovUnscaled)
})
