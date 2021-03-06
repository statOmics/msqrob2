% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accessors.R
\name{getModel,StatModel-method}
\alias{getModel,StatModel-method}
\alias{statModelAccessors}
\alias{getCoef}
\alias{getDF}
\alias{getDfPosterior}
\alias{getFitMethod}
\alias{getModel}
\alias{getSigma}
\alias{getSigmaPosterior}
\alias{getVar}
\alias{getVarPosterior}
\alias{getVcovUnscaled}
\alias{getFitMethod,StatModel-method}
\alias{getCoef,StatModel-method}
\alias{getDfPosterior,StatModel-method}
\alias{getVarPosterior,StatModel-method}
\alias{getSigmaPosterior,StatModel-method}
\alias{getDF,StatModel-method}
\alias{getVar,StatModel-method}
\alias{getSigma,StatModel-method}
\alias{getVcovUnscaled,StatModel-method}
\title{Accessor functions for StatModel class}
\usage{
\S4method{getModel}{StatModel}(object)

\S4method{getFitMethod}{StatModel}(object)

\S4method{getCoef}{StatModel}(object)

\S4method{getDfPosterior}{StatModel}(object)

\S4method{getVarPosterior}{StatModel}(object)

\S4method{getSigmaPosterior}{StatModel}(object)

\S4method{getDF}{StatModel}(object)

\S4method{getVar}{StatModel}(object)

\S4method{getSigma}{StatModel}(object)

\S4method{getVcovUnscaled}{StatModel}(object)
}
\arguments{
\item{object}{\code{StatModel} object}
}
\value{
The requested parameter of the StatModel object
}
\description{
Accessor functions for StatModel class
\describe{
\item{getModel(object)}{to get model}
\item{getFitMethod(object)}{to get the parameter estimation method}
\item{getCoef(object)}{to get the parameter estimates of the mean model}
\item{getDF(object)}{to get the residual degrees of freedom of the model}
\item{getVar(object)}{to get the residual variance of the model}
\item{getSigma(object)}{to get the residual standard deviation of the model}
\item{getDfPosterior(object)}{to get the degrees of freedom of
the empirical Bayes variance estimator}
\item{getVarPosterior(object)}{to get the empirical Bayes variance}
\item{getSigmaPosterior(object)}{to get the empirical Bayes standard deviation}
\item{getVcovUnscaled(object)}{to get the unscaled variance covariance matrix
of the model parameters}
}
}
\examples{
data(pe)

# Aggregate peptide intensities in protein expression values
pe <- aggregateFeatures(pe, i = "peptide", fcol = "Proteins", name = "protein")

# Fit msqrob model
pe <- msqrob(pe, i = "protein", formula = ~condition)
getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])
getModel(rowData(pe[["protein"]])$msqrobModels[[1]])
getFitMethod(rowData(pe[["protein"]])$msqrobModels[[1]])
# Similar for the remaining accessors
}
