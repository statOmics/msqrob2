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
{
vcovTmp<-getVcovUnscaled(object)*object@varPosterior
if (mean(rownames(vcovTmp)==rownames(L))==1)
return(t(L)%*%vcovTmp%*%L)
}
return(out)
}
)
