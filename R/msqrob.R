msqrob <- function(pe,assay,formula,nRobIter=1)
{
yAll<-assay(pe[[assay]])
ngenes<-nrow(yAll)
models<-apply(yAll,1,function(y,design)
{
        obs <- is.finite(y)
        if (sum(obs) > 0) {
            X <- design[obs, , drop = FALSE]
            y <- y[obs]
            mod<-lm.fit(X, y)
            modelType<-"ols"
            if (nRobIter>0)
            {
              for (kk in 1:nRobIter)
              {
                  w<-MASS::psi.huber(mod$res/mad(mod$res))
                  w[is.na(w)]<-0
                  mod<-lm.wfit(X, y,w=w)
                  if (mean(w)==1) break
              }
              modelType="robust"
            }
              out <- msqrobModel(modelType=modelType,model=mod)
            }
            else out<- msqrobModel(modelType="fitError",model=NULL)
            return(out)
},design=model.matrix(formula,colData(pe)))
hlp<-limma::squeezeVar(var=sapply(models,getVar),df=sapply(models,getDF))
out<-DataFrame(var.post=hlp$var.post,df.post=hlp$df.prior+sapply(models,getDF))
out$models<-models
return(out)
}
