#' Parameter estimates, standard errors and statistical inference on differential
#' expression
#'
#' @description Summary table of the estimates for differential expression of features
#'
#' @rdname hypothesisTest
#'
#' @aliases hypothesisTest hypothesisTest,SummarizedExperiment-method hypothesisTest,Features-method msqrobAggregate,Features-method
#'
#' @author Lieven Clement
#'
#' @examples
#'
#' # Load example data
#' # The data are a Feature object containing
#' # a SummarizedExperiment named "peptide" with MaxQuant peptide intensities
#' # The data are a subset of spike-in the human-ecoli study
#' # The variable condition in the colData of the Feature object
#' # contains information on the spike in condition a-e (from low to high)
#' data(pe)
#'
#' # Aggregate peptide intensities in protein expression values
#' pe<-aggregateFeatures(pe,i="peptide",fcol="Proteins",name="protein")
#'
#' # Fit msqrob model
#' pe <- msqrob(pe,i="protein",formula=~condition)
#'
#' # Define contrast
#' getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])
#' # Assess log2 fold change between condition c and condition b:
#' L <- makeContrast("conditionc - conditionb=0",c("conditionb","conditionc"))
#'
#' #example SummarizedExperiment instance
#' se <- pe[["protein"]]
#' se <- hypothesisTest(se,L)
#' head(rowData(se)$"conditionc - conditionb",10)
#' #Volcano plot
#' plot(-log10(pval)~logFC,rowData(se)$"conditionc - conditionb",col=(adjPval<0.05)+1)
#'
#' # Example for Features instance
#' # Assess log2 fold change between condition b and condition a (reference class),
#' # condition c and condition a, and, condition c and condition b.
#' L <-  makeContrast(c("conditionb=0","conditionc=0","conditionc - conditionb=0"),c("conditionb","conditionc"))
#' pe <- hypothesisTest(pe,i="protein",L)
#' head(rowData(pe[["protein"]])$"conditionb",10)
#' #Volcano plots
#' par(mfrow=c(1,3))
#' plot(-log10(pval)~logFC,rowData(pe[["protein"]])$"conditionb",col=(adjPval<0.05)+1,main="log2 FC b-a")
#' plot(-log10(pval)~logFC,rowData(pe[["protein"]])$"conditionc",col=(adjPval<0.05)+1,main="log2 FC c-a")
#' plot(-log10(pval)~logFC,rowData(pe[["protein"]])$"conditionc - conditionb",col=(adjPval<0.05)+1,main="log2 FC c-b")
#'
#' @param object `SummarizedExperiment` or `Features` instance
#' @param contrast `numeric` matrix specifying one or more contrasts of
#'        the linear model coefficients to be tested equal to zero. If multiple
#'        contrasts are given (multiple columns) then results will be returned for
#'        each contrast. The rownames of the matrix should be equal to the names
#'        of parameters of the model that are involved in the contrast.
#'        The column names of the matrix will be used to construct names to store
#'        the results in the rowData of the SummarizedExperiment or of the assay of
#'        the Features object. The contrast matrix can be made using the `makeContrast`
#'        function.
#' @param adjust.method `character` specifying the method to adjust
#'        the p-values for multiple testing.
#'        Options, in increasing conservatism, include ‘"none"’,
#'        ‘"BH"’, ‘"BY"’ and ‘"holm"’.  See ‘p.adjust’ for the complete
#'        list of options. Default is "BH" the Benjamini-Hochberg method
#'        to controle the False Discovery Rate (FDR).
#' @param modelColumn `character` to indicate the variable name that was used
#'        to store the msqrob models in the rowData of the SummarizedExperiment
#'        instance or of the assay of the Features instance. Default is "msqrobModels".
#' @param resultsColumnNamePrefix `character` to indicate the the prefix for the
#'        variable name that will be used to store test results in the rowData of
#'        the SummarizedExperiment instance or of the assay of the Features instance.
#'        Default is "" so that the variable name with the results will be
#'        the column name of the column in the contrast matrix L. If L is a matrix
#'        with multiple columns, multiple results columns will be made, one for each
#'        contrast. If L has no column names and if resultsColumnNamePrefix="" the
#'        results will be stored in the column with name msqrobResults.
#'
#' @export

setMethod("hypothesisTest","SummarizedExperiment",
          function(object,
                   contrast,
                   adjust.method="BH",
                   modelColumn="msqrobModels",
                   resultsColumnNamePrefix="",
                   overwrite=FALSE){
            if(!(modelColumn %in% colnames(rowData(object)))) stop(paste0("There is no column named \'", modelColumn,"\' with stored models of an msqrob fit in the rowData of the SummarizedExperiment object"))
            if(is.null(colnames(contrast)) & resultsColumnNamePrefix=="") resultsColumnNamePrefix<-"msqrobResults"
            if(is.null(colnames(contrast)) & ncol(contrast)>1) colnames(contrast) <- 1:ncol(contrast)
            if((sum(paste0(resultsColumnNamePrefix,colnames(contrast)) %in% colnames(rowData(object)))>0)&!overwrite) stop(paste0("There is/are already column(s) named \'", paste(paste0(resultsColumnNamePrefix,colnames(contrast)),collapse="\' or \'"),"\' in the rowData of the SummarizedExperiment object, set the argument overwrite=TRUE to replace the column(s) with the new results or use another name for the argument resultsColumnNamePrefix"))
            for (j in 1:ncol(contrast))
            {
                  contrHlp<-contrast[,j]
                  names(contrHlp)<-rownames(contrast)
                  rowData(object)[[paste0(resultsColumnNamePrefix,colnames(contrast)[j])]]<-topFeatures(rowData(object)[,modelColumn],contrast=contrHlp,adjust.method=adjust.method,sort=FALSE,alpha=1)
            }
            return(object)
            })

#' @param i `character` or `integer` to specify the element of the `Features` that
#'        contains the log expression intensities that will be modelled.
#'
#' @return A SummarizedExperiment or a `Features` instance augmented with the test
#'         results.
#'
#' @export
#' @rdname hypothesisTest

setMethod("hypothesisTest","Features",
          function(object,
                   i,
                   contrast,
                   adjust.method="BH",
                   modelColumn="msqrobModels",
                   resultsColumnNamePrefix="",
                   overwrite=FALSE){
            if (is.null(object[[i]])) stop(paste0("Features object does not contain an assay with the name ",i))
            if(!(modelColumn %in% colnames(rowData(object[[i]])))) stop(paste0("There is no column named \'", modelColumn,"\' with stored models of an msqrob fit in the rowData of assay ",i,"of the Features object."))
            if(is.null(colnames(contrast)) & resultsColumnNamePrefix=="") resultsColumnNamePrefix<-"msqrobResults"
            if(is.null(colnames(contrast)) & ncol(contrast)>1) colnames(contrast) <- 1:ncol(contrast)
            if((sum(paste0(resultsColumnNamePrefix,colnames(contrast)) %in% colnames(rowData(object[[i]])))>0)&!overwrite) stop(paste0("There is/are already column(s) named \'", paste(paste0(resultsColumnNamePrefix,colnames(contrast)),collapse="\' or \'"),"\' in the rowData of assay ",i," of the Features object, set the argument overwrite=TRUE to replace the column(s) with the new results or use another name for the argument resultsColumnNamePrefix"))
            for (j in 1:ncol(contrast))
            {
                  contrHlp<-contrast[,j]
                  names(contrHlp)<-rownames(contrast)
                  rowData(object[[i]])[[paste0(resultsColumnNamePrefix,colnames(contrast)[j])]]<-topFeatures(rowData(object[[i]])[,modelColumn],contrast=contrHlp,adjust.method=adjust.method,sort=FALSE,alpha=1)
            }
            return(object)
            })

#setMethod("Hurdle", "Features"
#          function(object,
#                   i,
#                   j,
#                   contrast,
#                   adjustMethod="BH",
#                   modelColumni="msqrobModels",
#                   modelColumnj="msqrobModels",
#                   resultsColumnNamePrefix=c("hurdleComp1","hurdleComp2","hurdle"))
#                   ){
#    object <- hypothesisTest(object,i,contrast,adjustMethod="BH",modelColumn=modelColumni,resultsColumnNamePrefix=resultsColumnNamePrefix[1])
#    object <- hypothesisTest(object,j,contrast,adjustMethod="BH",modelColumn=modelColumnj,resultsColumnNamePrefix=resultsColumnNamePrefix[1])
#    object[[i]]
#
#}
