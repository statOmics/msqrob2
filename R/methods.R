print.msqrobModel<-function(x,...)
{
cat("msqrobModel",x@modelType,"\n")

if (x@modelType=="ols"|x@modelType=="robust")
{
cat("Coefficients:\n")
cat(names(x@model$coefficients),"\n")
cat(x@model$coefficients,"\n")
}
}

setMethod("show","msqrobModel",function(object) print.msqrobModel(object))
