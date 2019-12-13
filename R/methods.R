setMethod("vcovUnscaled","StatModel",function(object)
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

setMethod("getContrast","StatModel",function(object,L)
{
coefs<-getCoef(object)
if(is.na(mean(names(coefs)==rownames(L))==1)) return(NA)
if (mean(names(coefs)==rownames(L))==1)
return(t(L)%*%coefs)
})

setMethod("varContrast","StatModel",function(object,L)
{
out<-matrix(NA,ncol(L),ncol(L))
rownames(out)<-colnames(out)<-colnames(L)
if (object@type!="fitError")
if (object@params$rank!=0)
{
vcovTmp<-vcovUnscaled(object)*object@varPosterior
if (mean(rownames(vcovTmp)==rownames(L))==1)
return(t(L)%*%vcovTmp%*%L)
}
return(out)
}
)

setMethod("getVar","StatModel",function(object)
{
mod<-getModel(object)
if(object@type=="fitError") return(NA)
if(object@type=="rlm")
{
return(sum(mod$w*mod$resid^2)/(sum(mod$w)-mod$rank))
} else if (is.null(mod$weights)) {
      return(sum(mod$residuals^2)/mod$df.residual)
        }
        else {
        return(sum(mod$weights*mod$residuals^2)/mod$df.residual)
        }
})
