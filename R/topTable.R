#' Toplist of DE proteins, peptides or features
#'
#' @description Summary table of the n most differentially expressed Features
#'
#' @param models A list with elements of the class StatModel that are estimated using the \code{\link{msqrob}} function
#' @param contrast A matrix with contrast. It is used to assess the null hypothesis that a linear combinations of the model parameters equals zero. The matrix needs to have an equal amount of rows as the number of parameters in the StatModel. The rownames should also match with the names of the model parameters/columns of the design matrix in the StatModel objects.
#' @examples TODO
#' @return A dataframe with log2 fold changes (logFC), standard errors (se), degrees of freedom of the test (df), t-test statistic (t), p-values (pval) and adjusted pvalues (adjPval) using the Benjamini-Hochberg method implemented in the p.adjust function of the stats package.
#' @rdname topTable
#' @author Lieven Clement
#' @export

topTable<-function(models,contrast)
{
logFC<-sapply(models,getContrast,L=contrast)
se<-sqrt(sapply(models,varContrast,L=contrast))
df<-sapply(models,getDfPosterior)
t<-logFC/se
pval<-pt(-abs(t),df)*2
adjPval<-p.adjust(pval,method="fdr")
return(data.frame(logFC,se,df,t,pval,adjPval))
}
