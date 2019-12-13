msqrob <- function(pe,assay,formula,robust=TRUE)
{
yAll<-assay(pe[[assay]])
ngenes<-nrow(yAll)
models<-apply(yAll,1,function(y,design)
{
        obs <- is.finite(y)
        modelType="fitError"
        model=list()
        if (sum(obs) > 0) {
            X <- design[obs, , drop = FALSE]
            y <- y[obs]
            if(robust)
            {
              model<-try(MASS::rlm(X, y,method="M"),silent=TRUE)
              if (class(model)[1]=="try-error") model=list() else
              {
                modelType<-"rlm"
                class(model)<-"list"
                }
            } else
            {
              model <- lm.fit(X, y)
              modelType <- "ols"
            }
            }
            return(msqrobModel(modelType=modelType,model=model))
},design=model.matrix(formula,colData(pe)))
hlp<-limma::squeezeVar(var=sapply(models,getVar),df=sapply(models,getDF))
for (i in 1:length(models))
{
  models[[i]]@varPosterior<-hlp$var.post[i]
  models[[i]]@dfPosterior<-hlp$df.prior+getDF(models[[i]])
}
return(models)
}
