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

setMethod("vcovUnscaled","msqrobModel",function(object)
{
    mod<-getModel(object)
    p1 <- 1L:mod$rank
    p <- length(mod$coef)
    out<-matrix(NA,p,p)
    out[p1,p1]<-chol2inv( mod$qr$qr[p1, p1, drop = FALSE])
    colnames(out)<-rownames(out)<-names(mod$coefficients)
    return(out)
}
)

setMethod("getContrast","msqrobModel",function(object,L)
{
coefs<-getCoef(object)
if(is.na(mean(names(coefs)==rownames(L))==1)) return(NA)
if (mean(names(coefs)==rownames(L))==1)
return(t(L)%*%coefs)
})

setMethod("varContrastUnscaled","msqrobModel",function(object,L)
{
if (object@model$rank==0) return(matrix(NA,1,1))
vcovTmp<-vcovUnscaled(object)
if (mean(rownames(vcovTmp)==rownames(L))==1)
return(t(L)%*%vcovTmp%*%L)
}
)

setMethod("getVar","msqrobModel",function(object)
{
mod<-getModel(object)
if(object@modelType=="fitError") return(NA)
if (is.null(mod$weights)) {
      return(sum(mod$residuals^2)/mod$df.residual)
        }
        else {
        return(sum(mod$weights*mod$residuals^2)/mod$df.residual)
        }
})
