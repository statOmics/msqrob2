##########################################################################################################
# Median Normalisation Factors
##########################################################################################################

#' Helper function to calculate sample-specific normalization factors on the log2 scale using 
#' conventional median normalisation
#'
#' @param mat A \code{matrix} object.
#' @param na.rm Logical; should missing values be removed? Default is \code{TRUE} as missing values typically occur in proteomics data.
#'
#' @return A numeric vector of log-scale normalization factors, one per sample (column).
#'
#' @details
#' This implementation assumes that the assay values are already on the log scale.
#' The normalization factors are computed on the log scale. 
#' @importFrom stats median


.computeNfLogMedian <- function(mat, na.rm = TRUE) {
  # 1. Calculates the sample medians of the intensity assay of summarised experiment `i` in qfeatures object `qf` 
  # 2. Defines the log norm factor by subtracting the median of the sample medians from each sample median,  so as to center the intensities of all samples around the median intensity in the experiment. 
  message("This function aims to calculate norm factors on a log scale, the input data are assumed to be on the log-scale!")
  nf_log <- mat |>  
    colMedians(na.rm = na.rm) #1.
  nf_log <- nf_log - median(nf_log) #2. 
  return(nf_log)
  }

#' Methods to computes sample-specific normalization factors on the log scale 
#' using conventional median summarisation.
#'
#' @aliases nfLogMedian nfLogMedian,SummarizedExperiment-method nfLogMedian,QFeatures-method nfLogMedian,matrix-method
#'
#' @param object A \code{matrix}, \code{SummarizedExperiment} or \code{QFeatures} object.
#' @param na.rm Logical; should missing values be removed? Default is \code{TRUE} as missing values typically occur in proteomics data.
#'
#' @return A numeric vector of log2-scale normalization factors, one per sample (column).
#'
#' @details
#' This implementation assumes that the assay values are already on the log scale.
#' The normalization factors are computed on the log scale. 
#'
#' @examples
#' # Load example data
#' # The data are a Feature object with containing
#' # a SummarizedExperiment named "peptide" with MaxQuant peptide intensities
#' # The data are a subset of spike-in the human-ecoli study
#' # The variable condition in the colData of the Feature object
#' # contains information on the spike in condition a-e (from low to high)
#' data(pe)
#' 
#' ###############################
#' ### Example on QFeatures object
#' ###############################
#' 
#' data(pe)
#'
#' # Calculate log2 norm factor
#' nf_log <- nfLogMedian(pe, i="peptide")
#' nf_log
#' 
#' # Normalise peptide level data
#' pe <- sweep(pe, 
#'     MARGIN = 2, 
#'     STATS = nf_log, 
#'     i = "peptide", 
#'     name = "peptide_norm")
#'     
#' # Evaluate normalisation
#' par(mfrow=c(1,2))
#' boxplot(assay(pe,"peptide"))
#' boxplot(assay(pe,"peptide_norm"))
#' 
#' ###############################
#' ### Example on SummarizedExperiment object
#' ###############################
#' 
#' data(pe)
#'
#' # Extract a summarised experiment from QFeatures object pe
#' se <- getWithColData(pe, i="peptide")
#' 
#' # Calculate log2 norm factor
#' nf_log <- nfLogMedian(se, i=1)
#' nf_log
#' 
#' # Normalise peptide level data and store it as a new assay in
#' # the SummarizedExperiment object
#' assays(se)[["peptide_norm"]] <- sweep(assay(se), 
#'     MARGIN = 2, 
#'     STATS = nf_log) 
#'     
#' # Also give first assay a name     
#' assayNames(se)[1] <- "peptide"
#' 
#' # Evaluate normalisation
#' par(mfrow=c(1,2))
#' boxplot(assay(se,"peptide"))
#' boxplot(assay(se,"peptide_norm"))
#' 
#' ###############################
#' ### Example on matrix object
#' ###############################
#'
#' data(pe)
#'
#' # Extract log2 transformed intensity  matrix from QFeatures object pe
#' mat <- assay(pe,"peptide")
#' 
#' # Calculate Norm factors
#' nf <- nfLogMedian(mat)
#' nf
#' 
#' # Normalise peptide level data
#' matnorm <- sweep(mat, 
#'     MARGIN = 2, 
#'     STATS = nf_log) 
#' # Evaluate normalisation
#' par(mfrow=c(1,2))
#' boxplot(mat)
#' boxplot(matnorm)
#'
#' @importFrom matrixStats colMedians
#' @importFrom stats median
#' @export
#' @rdname nfLogMedian
setMethod("nfLogMedian", signature(object = "matrix"),
          function(object, na.rm = TRUE) {
            .computeNfLogMedian(mat = object, na.rm=na.rm)
          })

#' @param i An integer or character specifying which assay to use, only needed when object is \code{SummarizedExperiment} or \code{QFeatures}
#' @export
#' @rdname nfLogMedian
setMethod("nfLogMedian", signature(object = "SummarizedExperiment"),
          function(object, i, na.rm = TRUE) {
            if (missing(i)) stop("No assay provided, please define argument 'i'")
            if (class(try(object[[i]], silent=TRUE)) %in% c("try-error","NULL")) stop("Object does not contain an assay with the name ", i)
            .computeNfLogMedian(mat = SummarizedExperiment::assay(object, i), na.rm=na.rm)
          })

#' @param i An integer or character specifying which assay to use, only needed when object is \code{SummarizedExperiment} or \code{QFeatures}
#' @export
#' @rdname nfLogMedian
setMethod("nfLogMedian", signature(object = "QFeatures"),
          function(object, i, na.rm = TRUE) {
            if (missing(i)) stop("No assay provided, please define argument 'i'")
            if (class(try(object[[i]], silent=TRUE)) %in% c("try-error","NULL")) stop("Object does not contain an assay with the name ", i)
            .computeNfLogMedian(mat = SummarizedExperiment::assay(object, i), na.rm=na.rm)
          })

##########################################################################################################
# Median-of-Ratios Normalisation Factors
##########################################################################################################

#' Helper function to calculate sample-specific normalization factors on the log2 scale using a
#' median-of-ratios approach similar to that used in DESeq2 for bulk RNA-seq data.
#'
#' The method proceeds as follows:
#' \enumerate{
#'   \item A pseudo-reference sample is constructed as the row-wise mean of the
#'   log2 intensities (equivalent to the log2-transformed geometric mean).
#'   \item For each sample, log2 ratios relative to the pseudo-reference are computed.
#'   \item The normalization factor for each sample is obtained as the median of
#'   these log2 ratios (column-wise median).
#' }
#' @param mat A \code{matrix} object.
#' @param na.rm Logical; should missing values be removed? Default is \code{TRUE} as missing values typically occur in proteomics data.
#'
#' @return A numeric vector of log2-scale normalization factors, one per sample (column).
#'
#' @details
#' This implementation assumes that the assay values are already on the log scale.
#' The normalization factors are computed on the log scale. 

.computeNfLogMedianOfRatios <- function(mat, na.rm = TRUE) {
  #1. Calculate reference sample using the feature mean 
  #   as the data are on log-scale this is the log-transformed geometric mean
  #2. Calculate logFC w.r.t. the reference sample. 
  #   By default FUN argument in sweep is "-" (subtract)
  #   Margin = 1 --> adopt FUN on rows
  #3. Calculate median logFC per sample, which is the final log scale norm factor
  message("This function aims to calculate norm factors on a log scale, the input data are assumed to be on the log-scale!")
  pseudoref <- rowMeans(mat, na.rm = na.rm) #1. 
  sweep(mat, MARGIN = 1, pseudoref) |>      #2.
    matrixStats::colMedians(na.rm = na.rm)  #3.
}


#' Methods to calculate log-scale normalization factors using the median-of-ratios method
#'
#' @description Computes sample-specific normalization factors on the log2 scale using a
#' median-of-ratios approach similar to that used in DESeq2 for bulk RNA-seq data.
#'
#' The method proceeds as follows:
#' \enumerate{
#'   \item A pseudo-reference sample is constructed as the row-wise mean of the
#'   log2 intensities (equivalent to the log2-transformed geometric mean).
#'   \item For each sample, log2 ratios relative to the pseudo-reference are computed.
#'   \item The normalization factor for each sample is obtained as the median of
#'   these log2 ratios (column-wise median).
#' }
#'
#' @aliases nfLogMedianOfRatios nfLogMedianOfRatios,SummarizedExperiment-method nfLogMedianOfRatios,QFeatures-method nfLogMedianOfRatios,matrix-method
#' @param object A \code{matrix}, \code{SummarizedExperiment} or \code{QFeatures} object.
#' @param na.rm Logical; should missing values be removed? Default is \code{TRUE} as missing values typically occur in proteomics data.
#'
#' @return A numeric vector of log2-scale normalization factors, one per sample (column).
#'
#' @details
#' This implementation assumes that the assay values are already on the log scale.
#' The normalization factors are computed on the log scale. 
#'
#' @examples
#' # Load example data
#' # The data are a Feature object with containing
#' # a SummarizedExperiment named "peptide" with MaxQuant peptide intensities
#' # The data are a subset of spike-in the human-ecoli study
#' # The variable condition in the colData of the Feature object
#' # contains information on the spike in condition a-e (from low to high)
#' data(pe)
#' 
#' ###############################
#' ### Example on QFeatures object
#' ###############################
#' 
#' data(pe)
#'
#' # Calculate log2 norm factor
#' nf_log <- nfLogMedianOfRatios(pe, i="peptide")
#' nf_log
#' 
#' # Normalise peptide level data
#' pe <- sweep(pe, 
#'     MARGIN = 2, 
#'     STATS = nf_log, 
#'     i = "peptide", 
#'     name = "peptide_norm")
#'     
#' # Evaluate normalisation
#' par(mfrow=c(1,2))
#' boxplot(assay(pe,"peptide"))
#' boxplot(assay(pe,"peptide_norm"))
#' 
#' ###############################
#' ### Example on SummarizedExperiment object
#' ###############################
#' 
#' data(pe)
#'
#' # Extract a summarised experiment from QFeatures object pe
#' se <- getWithColData(pe, i="peptide")
#' 
#' # Calculate log2 norm factor
#' nf_log <- nfLogMedianOfRatios(se, i=1)
#' nf_log
#' 
#' # Normalise peptide level data and store it as a new assay in
#' # the SummarizedExperiment object
#' assays(se)[["peptide_norm"]] <- sweep(assay(se), 
#'     MARGIN = 2, 
#'     STATS = nf_log) 
#'     
#' # Also give first assay a name     
#' assayNames(se)[1] <- "peptide"
#' 
#' # Evaluate normalisation
#' par(mfrow=c(1,2))
#' boxplot(assay(se,"peptide"))
#' boxplot(assay(se,"peptide_norm"))
#' 
#' ###############################
#' ### Example on matrix object
#' ###############################
#' 
#' data(pe)
#'
#' mat <- assay(pe,"peptide")
#' 
#' nf <- nfLogMedianOfRatios(mat)
#' nf
#' 
#' # Normalise peptide level data
#' matnorm <- sweep(mat, 
#'     MARGIN = 2, 
#'     STATS = nf_log) 
#' # Evaluate normalisation
#' par(mfrow=c(1,2))
#' boxplot(mat)
#' boxplot(matnorm)
#'
#' @references
#' Love, M.I., Huber, W., Anders, S. (2014).
#' Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.
#' Genome Biology, 15(12), 550.
#'
#' @importFrom matrixStats colMedians
#' @export
#' @rdname nfLogMedianOfRatios
setMethod("nfLogMedianOfRatios", signature(object = "matrix"),
          function(object, na.rm = TRUE) {
            .computeNfLogMedianOfRatios(mat = object, na.rm=na.rm)
          })

#' @param i An integer or character specifying which assay to use, only needed when object is \code{SummarizedExperiment} or \code{QFeatures}
#' @export
#' @rdname nfLogMedianOfRatios
setMethod("nfLogMedianOfRatios", signature(object = "SummarizedExperiment"),
          function(object, i, na.rm = TRUE) {
            if (missing(i)) stop("No assay provided, please define argument 'i'")
            if (class(try(object[[i]], silent=TRUE)) %in% c("try-error","NULL")) stop("Object does not contain an assay with the name ", i)
            .computeNfLogMedianOfRatios(mat = SummarizedExperiment::assay(object, i), na.rm=na.rm)
          })

#' @param i An integer or character specifying which assay to use, only needed when object is \code{SummarizedExperiment} or \code{QFeatures}
#' @export
#' @rdname nfLogMedianOfRatios
setMethod("nfLogMedianOfRatios", signature(object = "QFeatures"),
          function(object, i, na.rm = TRUE) {
            if (missing(i)) stop("No assay provided, please define argument 'i'")
            if (class(try(object[[i]], silent=TRUE)) %in% c("try-error","NULL")) stop("Object does not contain an assay with the name ", i)
            .computeNfLogMedianOfRatios(mat = SummarizedExperiment::assay(object, i), na.rm=na.rm)
          })
