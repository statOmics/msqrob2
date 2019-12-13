#' Function to fit msqrob models
#'
#' @description Parameter estimation of msqrob models.
#'
#' @param pe Proteomics Experiment Object, is a Feature set that contains the assay with the quantified MS intensities.
#' @param assay integer or string referring to the assay of the Feature object containing the intensities for modelling, e.g. peptide or protein.
#' @param formula Model formula. The model is build based on the covariates in the colData of the pe object.
#' @param robust Flag to indicate if robust regression is conducted to account for outliers. If Flase an OLS fit is performed.
#' @param maxitRob Maximum iterations in the IRWLS algorithm used in the M-estimation step of the robust regression.
#' @examples TODO
#' @return A list of objects of the StatModel class
#' @rdname msqrob
#' @author Lieven Clement
#' @export
msqrob <- function(pe,assay,formula,robust=TRUE,maxitRob=5)
{
yAll<-assay(pe[[assay]])
ngenes<-nrow(yAll)
models<-apply(yAll,1,function(y,design)
{
        obs <- is.finite(y)
        type="fitError"
        model=list()
        if (sum(obs) > 0) {
            X <- design[obs, , drop = FALSE]
            y <- y[obs]
            if(robust)
            {
              model<-try(MASS::rlm(X, y,method="M",maxit=maxitRob),silent=TRUE)
              if (class(model)[1]=="try-error") model=list() else
              {
                type<-"rlm"
                class(model)<-"list"
                }
            } else
            {
              model <- lm.fit(X, y)
              type <- "lm"
            }
            }
            return(StatModel(type=type,params=model,varPosterior=as.numeric(NA),dfPosterior=as.numeric(NA)))
},design=model.matrix(formula,colData(pe)))
hlp<-limma::squeezeVar(var=sapply(models,getVar),df=sapply(models,getDF))
for (i in 1:length(models))
{
  models[[i]]@varPosterior<-as.numeric(hlp$var.post[i])
  models[[i]]@dfPosterior<-as.numeric(hlp$df.prior+getDF(models[[i]]))
}
return(models)
}
