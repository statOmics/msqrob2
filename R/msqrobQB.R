
#' Function to fit msqrob models to peptide counts using glm
#'
#' @description Low-level function for parameter estimation with msqrob
#'              by modeling peptide counts using quasibinomial glm
#'
#' @param object `SummarizedExperiment` or `Features` instance
#'
#' @param formula Model formula. The model is built based on the
#'        covariates in the data object.
#'
#' @param modelColumnName `character` to indicate the variable name that is used
#'        to store the msqrob models in the rowData of the SummarizedExperiment
#'        instance or of the assay of the Features instance. Default is "msqrobModels".
#'
#' @param overwrite `boolean(1)` to indicate if the column in the rowData has to
#'        be overwritten if the modelColumnName already exists. Default is FALSE.
#'
#' @param nFeat `character(1)` to indicate the name of the variable in the rowData
#'        of object that contains the maximum number of features (maximum count).
#'        E.g. when peptides are summarized to proteins by counting how many peptides
#'        that were observed for a protein in a sample, nFeat will indicate the maximum
#'        number of peptides that could be observed for a protein. The default is "nPep"
#'
#' @param priorCount A 'numeric(1)', which is a prior count to be added to the observations to shrink
#'          the estimated log-fold-changes towards zero. Default is 0.1.
#'
#' @param binomialBound: logical, if ‘TRUE’ then the quasibinomial variance estimator will
#'                be never smaller than 1 (no underdispersion). Default is TRUE.
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
#' # Aggregate by counting how many peptide we observe for each protein
#' pe<-aggregateFeatures(pe,i="peptide",fcol="Proteins",name="proteinCount",fun=nObsPep)
#'
#' # Calculate the total number of peptides for each protein
#' rowData(pe[["proteinCount"]])$nPep<- rowData(pe[["peptide"]]) %>%
#'        data.frame(.$Proteins)  %>%
#'        group_by(Proteins) %>%
#'        summarise(nPep = length(Proteins)) %>%
#'        column_to_rownames("Proteins") %>%
#'        .$nPep
#'
#' # Fit MSqrob model to peptide counts using a quasi-binomial model
#' # For summarized SummarizedExperiment
#' se <- pe[["proteinCount"]]
#' se
#' colData(se) <- colData(pe)
#' se <- msqrobQB(se,formula=~condition)
#' getCoef(rowData(se)$msqrobModels[[1]])
#'
#' # For features object
#' pe <- msqrobQB(pe,i="proteinCount",formula=~condition)
#'
#' @return SummarizedExperiment or Features instance
#'
#' @aliases msqrobQB msqrobQB,SummarizedExperiment-method msqrobQB,Features-method
#'
#' @author Lieven Clement
#'
#' @export

setMethod("msqrobQB","SummarizedExperiment",
          function(object,
                   formula,
                   modelColumnName="msqrobModels",
                   overwrite=FALSE,
                   nFeat="nPep",
                   priorCount=.1,
                   binomialBound=TRUE){

           if (ncol(colData(object))==0) stop("error: colData is empty")
           if((modelColumnName %in% colnames(rowData(object)))&!overwrite) stop(paste0("There is already a column named \'",
                                                                               modelColumnName,
                                                                               "\' in the rowData of the SummarizedExperiment object, set the argument overwrite=TRUE to replace the column with the new results or use another name for the argument modelColumnName to store the results as a novel column in the rowData of SummarizedExperiment object"))
           if(!(nFeat %in% colnames(rowData(object)))) stop(paste0("There is no column named\'",
                                                                   nFeat,
                                                                   "\'in the rowData of object"))
           rowData(object)[[modelColumnName]] <- msqrobGlm(assay(object),
                                                           rowData(object)[[nFeat]],
                                                           formula,
                                                           colData(object),
                                                           priorCount=priorCount,
                                                           binomialBound=binomialBound)
           return(object)
})

#' @param i `character` or `integer` to specify the element of the `Features` that
#'        contains the log expression intensities that will be modelled.
#' @export
#' @rdname msqrobQB

setMethod("msqrobQB","Features",
          function(object,
                   i,
                   formula,
                   modelColumnName="msqrobModels",
                   overwrite=FALSE,
                   nFeat="nPep",
                   priorCount=.1,
                   binomialBound=TRUE){
           if (is.null(object[[i]])) stop(paste0("Features object does not contain an assay with the name ",i))
           if((modelColumnName %in% colnames(rowData(object[[i]])))&!overwrite) stop(paste0("There is already a column named \'",
                                                                               modelColumnName,
                                                                               "\' in the rowData of assay \'",
                                                                               i,
                                                                               "'of object, set the argument overwrite=TRUE to replace the column with the new results or use another name for the argument modelColumnName to store the results as a novel column in the rowData of SummarizedExperiment object"))
           if(!(nFeat %in% colnames(rowData(object[[i]])))) stop(paste0("There is no column named\'",
                                                                   nFeat,
                                                                   "\'in the rowData of assay \'",
                                                                   i,
                                                                   "\' of object"))
           rowData(object[[i]])[[modelColumnName]] <- msqrobGlm(assay(object[[i]]),
                                                           rowData(object[[i]])[[nFeat]],
                                                           formula,
                                                           colData(object),
                                                           priorCount=priorCount,
                                                           binomialBound=binomialBound)
           return(object)
})
