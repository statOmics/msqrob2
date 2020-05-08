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
#' # The data are a Feature object with containing
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
#' #Define contrast
#' coef1<-getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])
#' L <- matrix(0,ncol=1,nrow=length(coef1))
#' rownames(L) <- names(coef1)
#' L[2,1]<-1
#'
#' #example SummarizedExperiment instance
#' se <- pe[["protein"]]
#' se <- hypothesisTest(se,L)
#' head(rowData(se)$msqrobTestResults,10)
#'
#' #example Features instance
#' pe <- hypothesisTest(pe,i="protein",L)
#' head(rowData(pe[["protein"]])$msqrobTestResults,10)
#' #Volcano plot
#' plot(-log10(pval)~logFC,rowData(pe[["protein"]])$msqrobTestResults,col=(adjPval<0.05)+1)
#' @param object `SummarizedExperiment` or `Features` instance
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
#' @param modelColumn `character` to indicate the variable name that is used
#'        to store the msqrob models in the rowData of the SummarizedExperiment
#'        instance or of the assay of the Features instance. Default is "msqrobModels".
#' @param resultsColumnName `character` to indicate the variable name that will be used
#'        to store test results in the rowData of the SummarizedExperiment instance
#'        or of the assay of the Features instance. Default is "msqrobTestResults".
#'
#' @export

setMethod("hypothesisTest","SummarizedExperiment",
          function(object,
                   contrast,
                   adjust.method="BH",
                   modelColumn="msqrobModels",
                   resultsColumnName="msqrobTestResults",
                   overwrite=FALSE){
            if(!(modelColumn %in% colnames(rowData(object)))) stop(paste0("There is no column named \'", modelColumn,"\' with stored models of an msqrob fit in the rowData of the SummarizedExperiment object"))
            if((resultsColumnName %in% colnames(rowData(object)))&!overwrite) stop(paste0("There is already a column named \'", resultsColumnName,"\' in the rowData of the SummarizedExperiment object, set the argument overwrite=TRUE to replace the column with the new results or use another name for the argument resultsColumnName to store the results as a novel column in the rowData of the SummarizedExperiment"))
            rowData(object)[[resultsColumnName]]<-topFeatures(rowData(object)[,modelColumn],contrast=contrast,adjust.method=adjust.method,sort=FALSE,alpha=1)
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
                   resultsColumnName="msqrobTestResults",
                   overwrite=FALSE){
              if (is.null(object[[i]])) stop(paste0("Features object does not contain an assay with the name ",i))
              if(!(modelColumn %in% colnames(rowData(object[[i]])))) stop(paste0("There is no column named \'", modelColumn,"\' with stored models of an msqrob fit in the rowData of the chosen assay of the Features object"))
              if((resultsColumnName %in% colnames(rowData(object[[i]])))&!overwrite) stop(paste0("There is already a column named \'", resultsColumnName,"\' in the rowData of the chosen assay of the Features object, set the argument overwrite=TRUE to replace the column with the new results or use another name for the argument resultsColumnName to store the results as a novel column in the rowData of assay"))
              rowData(object[[i]])[[resultsColumnName]]<-topFeatures(rowData(object[[i]])[,modelColumn],contrast=contrast,adjust.method=adjust.method,sort=FALSE,alpha=1)
              return(object)
              })
