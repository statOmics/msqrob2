#' Smallest unique protein groups
#'
#' @description For a given vector of protein group names, outputs the names of those protein groups for which none of its member proteins is present in a smaller protein group.
#' @param proteins A vector of characters or factors containing single proteins and/or protein groups (i.e. proteins separated by a separator symbol).
#' @param split The character string that is used to separate the indivudual protein names in each protein group.
#' @return A character vector containing the names of the protein groups for which none of its proteins is present in a smaller protein group.
#' @examples #TODO
#' @export
smallestUniqueGroups <- function(proteins,
                                 split=";"){
                        b <- strsplit(x=as.character(unique(proteins)),split=split,fixed=TRUE)

                        included <- vector()

                        j <- 1
                        while(length(b)!=0){
                            included <- c(included,sapply(b[sapply(b, length)==j], function(x) paste(x, collapse=split)))
                            a <- unlist(b[sapply(b, length)==j])
                            b <- b[sapply(b, length)>j]

                            if(length(b)!=0){
                                sel <- vector()
                                for(i in 1:length(b)){
                                    sel[i] <- !any(b[[i]] %in% a)
                                }
                            b <- b[sel]
                            j <- j+1
                            }
                        }

                        included <- unlist(included)
                        return(included)
                        }


#' Make contrast matrix
#'
#' @description  Construct the contrast matrix corresponding to specified contrasts
#'               of a set of parameters.
#'
#' @param contrasts: character vector specifying contrasts, i.e. the linear combination of the modelparameters that equals to zero.
#'
#' @param parameterNames: character vector specifying the model parameters that are involved in the contrasts, e.g if we model data of
#'        three conditions using a factor condition with three levels a, b and c then our model will have 3 mean parameters named (Intercept),
#'        conditionb and conditionc. Hence the log2 fold change between  b and a is conditionb. Under the null hypothesis the log2 fold change
#'        equals 0. Which is to be encoded as "conditionb=0". If we would like to test for log2 fold change between condition c and b we assess if
#'        the log2 fold change conditionc-conditionb equals 0, encoded as "condtionb-conditionc=0".
#'
#' @examples
#' makeContrast(c("conditionb=0"),parameterNames=c("(Intercept)","conditionb","conditionc"))
#' makeContrast(c("conditionc=b"),parameterNames=c("conditionB"))
#' makeContrast(c("conditionb=0","conditionc=0","conditionc-conditionb=0"),parameterNames=c("conditionb","conditionc"))
#'
#' @return A numeric contrast matrix with rownames that equal the model parameters that are involved in the contrasts
#'
#' @rdname makeContrast
#'
#' @export


makeContrast<-function(contrasts,parameterNames){
                  return(t(multcomp:::chrlinfct2matrix(contrasts,parameterNames)$K))
              }


#' Summary function for aggregateFeatures to count the number of observed peptides
#' per sample
#'
#' @description  Summarise function for aggregateFeatures to count the number of
#'               observed features.
#'
#' @param y Matrix of features
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
#' # Aggregate peptide intensities in protein expression values by counting number
#' # of peptides per protein
#' pe<-aggregateFeatures(pe,i="peptide",fcol="Proteins",name="proteinCount",fun=nObsPep)
#'
#' #Add total number of peptides per Protein in the rowData of the new assay
#' rowData(pe[["proteinCount"]])$nPep <- rowData(pe[["peptide"]]) %>%
#'        data.frame(.$Proteins)  %>%
#'        group_by(Proteins) %>%
#'        summarise(no_rows = length(Proteins)) %>%
#'        column_to_rownames("Proteins") %>%
#'        .$(nPep)
#'
#' @return A vector with counts
#'
#' @rdname nPep
#'
#' @export

nObsPep <- function(y) colSums(!is.na(y))
