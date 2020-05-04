setMethod("getContrast","StatModel",function(object,L)
{
coefs<-getCoef(object)
out<-matrix(rep(NA,ncol(L)))
rownames(out)<-colnames(L)
hlp<- try(t(L)%*%coefs[rownames(L)],silent=TRUE)
if (class(out)[1]!="try-error") out[]<-hlp
return(out)
})

setMethod("varContrast","StatModel",function(object,L)
{
out<-matrix(NA,ncol(L),ncol(L))
rownames(out)<-colnames(out)<-colnames(L)
vcovTmp<-getVcovUnscaled(object)*object@varPosterior
hlp<-try(t(L)%*%vcovTmp[rownames(L),rownames(L)]%*%L,silent=TRUE)
if (class(hlp)[1]!="try-error") out[]<-hlp
return(out)
}
)

setMethod("msqrob","SummarizedExperiment",
  function(object,formula,modelColumnName="msqrobModels",overwrite=FALSE,robust=TRUE,maxitRob=1,ridge=FALSE)
  {
   if (ncol(colData(object))==0) stop("error: colData is empty")
   if((modelColumnName %in% colnames(rowData(object)))&!overwrite) stop(paste0("There is already a column named \'", modelColumnName,"\' in the rowData of the SummarizedExperiment object, set the argument overwrite=TRUE to replace the column with the new results or use another name for the argument modelColumnName to store the results as a novel column in the rowData of SummarizedExperiment object"))
   if (!ridge)  rowData(object)[[modelColumnName]]<-msqrobLm(y=assay(object),formula=formula,data=colData(object),robust=robust,maxitRob=maxitRob) else
                rowData(object)[[modelColumnName]]<-msqrobLmer(y=assay(object),formula=formula,data=colData(object),robust=robust,maxitRob=maxitRob)
   return(object)
  }
)

setMethod("hypothesisTest","SummarizedExperiment",function(object,contrast,modelColumn="msqrobModels",contrastName="msqrobTestResults",overwrite=FALSE){
if(!(modelColumn %in% colnames(rowData(object)))) stop(paste0("There is no column named \'", modelColumn,"\' with stored models of an msqrob fit in the rowData of the SummarizedExperiment object"))
if((contrastName %in% colnames(rowData(object)))&!overwrite) stop(paste0("There is already a column named \'", contrastName,"\' in the rowData of the SummarizedExperiment object, set the argument overwrite=TRUE to replace the column with the new results or use another name for the argument contrastName to store the results as a novel column in the rowData of the SummarizedExperiment"))
rowData(object)[[contrastName]]<-topTable(rowData(object)[,modelColumn],contrast)
return(object)
})

setMethod("msqrob","Features",
  function(object,formula,assayName,modelColumnName="msqrobModels",overwrite=FALSE,robust=TRUE,maxitRob=1,ridge=FALSE)
  {
   if (is.null(object[[assayName]])) stop(paste0("Features object does not contain an assay with the name ",assayName))
   if((modelColumnName %in% colnames(rowData(object[[assayName]])))&!overwrite) stop(paste0("There is already a column named \'", modelColumnName,"\' in the rowData of the assay",assayName,"of the Features object, set the argument overwrite=TRUE to replace the column with the new results or use another name for the argument modelColumnName to store the results as a novel column in the rowData of assay of the Features object"))
   if (!ridge)  rowData(object[[assayName]])[[modelColumnName]]<-msqrobLm(y=assay(object[[assayName]]),formula=formula,data=colData(object),robust=robust,maxitRob=maxitRob) else
                rowData(object[[assayName]])[[modelColumnName]]<-msqrobLmer(y=assay(object[[assayName]]),formula=formula,data=colData(object),robust=robust,maxitRob=maxitRob)
   return(object)
  }
)

setMethod("hypothesisTest","Features",function(object,assayName,contrast,modelColumn="msqrobModels",contrastName="msqrobTestResults",overwrite=FALSE){
if (is.null(object[[assayName]])) stop(paste0("Features object does not contain an assay with the name ",assayName))
if(!(modelColumn %in% colnames(rowData(object[[assayName]])))) stop(paste0("There is no column named \'", modelColumn,"\' with stored models of an msqrob fit in the rowData of the chosen assay of the Features object"))
if((contrastName %in% colnames(rowData(object[[assayName]])))&!overwrite) stop(paste0("There is already a column named \'", contrastName,"\' in the rowData of the chosen assay of the Features object, set the argument overwrite=TRUE to replace the column with the new results or use another name for the argument contrastName to store the results as a novel column in the rowData of assay"))
rowData(object[[assayName]])[[contrastName]]<-topTable(rowData(object[[assayName]])[,modelColumn],contrast)
return(object)
})
