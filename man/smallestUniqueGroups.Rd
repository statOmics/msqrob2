% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msqrob-utils.R
\name{smallestUniqueGroups}
\alias{smallestUniqueGroups}
\title{Smallest unique protein groups}
\usage{
smallestUniqueGroups(proteins, split = ";")
}
\arguments{
\item{proteins}{A vector of characters or factors containing single proteins and/or protein groups (i.e. proteins separated by a separator symbol).}

\item{split}{The character string that is used to separate the indivudual protein names in each protein group.}
}
\value{
A character vector containing the names of the protein groups for which none of its proteins is present in a smaller protein group.
}
\description{
For a given vector of protein group names, outputs the names of those protein groups for which none of its member proteins is present in a smaller protein group.
}
\examples{
data(pe)
smallestUniqueGroups(rowData(pe[["peptide"]])$Proteins)
}
