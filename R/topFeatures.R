#' Toplist of DE proteins, peptides or features
#'
#' @description Summary table of the differentially expressed Features
#'
#' @param models A list with elements of the class `StatModel` that are
#'        estimated using the \code{\link{msqrob}} function
#' @param contrast `numeric` matrix specifying one or more contrasts of
#'        the linear model coefficients to be tested equal to zero.
#'        The rownames of the matrix should be equal to the names of
#'        parameters of the model.
#' @param adjust.method `character` specifying the method to adjust
#'        the p-values for multiple testing.
#'        Options, in increasing conservatism, include ‘"none"’,
#'        ‘"BH"’, ‘"BY"’ and ‘"holm"’.  See ‘p.adjust’ for the complete
#'        list of options. Default is "BH" the Benjamini-Hochberg method
#'        to controle the False Discovery Rate (FDR).
#' @param `boolean(1)` to indicate if the features have to be sorted according
#'        to statistical significance.
#' @param alpha `numeric` specifying the cutoff value for adjusted p-values.
#'        Only features with lower p-values are listed.
#'
#' @examples #TODO
#' @return A dataframe with log2 fold changes (logFC), standard errors (se),
#'         degrees of freedom of the test (df), t-test statistic (t),
#'         p-values (pval) and adjusted pvalues (adjPval) using the specified
#'         adjust.method in the p.adjust function of the stats package.
#' @rdname topTable
#' @author Lieven Clement
#' @export

topFeatures<-function(models,contrast,adjust.method="BH",sort=TRUE,alpha=1){
    logFC<-sapply(models,getContrast,L=contrast)
    se<-sqrt(sapply(models,varContrast,L=contrast))
    df<-sapply(models,getDfPosterior)
    t<-logFC/se
    pval<-pt(-abs(t),df)*2
    adjPval<-p.adjust(pval,method=adjust.method)
    out <-data.frame(logFC,se,df,t,pval,adjPval)
    if (sort) {
        if(alpha<1) {
            ids <- adjPval<alpha
            out<-na.exclude(out[ids,])
            }
        return(out[order(out$pval),])
        } else {
        return(out)
        }
}
