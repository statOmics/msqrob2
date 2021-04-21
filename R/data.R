#' Example data for 100 proteins
#' @description
#' Subset of peptides from 100 proteins from a quantitative mass spectrometry
#' based proteomics dataset (PRIDE identifier: PXD003881 Shen et al. (2018)).
#' E. Coli lysates were spiked at five different concentrations (3%, 4.5%, 6%,
#' 7.5% and 9%wt/wt) in a stable human background (4 repl. per treatment). The
#' twenty resulting samples were run on an Orbitrap Fusion mass spectrometer.
#' Raw data files were processed with MaxQuant (version 1.6.1.0, Cox and Mann
#' (2008)) using default search settings unless otherwise noted. Spectra were
#' searched against the UniProtKB/SwissProt human and E. Coli reference
#' proteome databases (07/06/2018), concatenated with the default Maxquant
#' contaminant database. Carbamidomethylation of Cystein was set as a fixed
#' modification, and oxidation of Methionine and acetylation of the protein
#' amino-terminus were allowed as variable modifications. In silico cleavage
#' was set to use trypsin/P, allowing two miscleavages. Match between runs
#' was also enabled using default settings. The resulting peptide-to-spectrum
#' matches (PSMs) were filtered by MaxQuant at 1% FDR.
#'
#' @format Feature set with an instance "peptide":
#'         \describe{
#'         \item{assay}{contains the raw peptide intensities}
#'         \item{rowData}{contains a variable "Proteins" with the protein accession and an variable ecoli to indicate if the protein is a spikin.}
#'         \item{colData}{contains a factor condition indicating the spike-in condition}
#'         }
#'
#' @rdname data
#'
#' @aliases data pe
#'
#' @examples
#' data(pe)
#' head(colData(pe))
#' head(rowData(pe))
#' head(assay(pe))
"pe"
