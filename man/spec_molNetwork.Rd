% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spec.R
\name{spec_molNetwork}
\alias{spec_molNetwork}
\title{Calculate molecular network from MS2 data}
\usage{
spec_molNetwork(x, methods = c("ndotproduct"), ...)
}
\arguments{
\item{x}{\code{Spectra} Spectra object with single spectrum per feature}

\item{methods}{\code{character} Name of function(s) used for comparison}

\item{...}{Parameters used by compareSpectra and sub functions}
}
\value{
\code{list} a named list with adjency matrices for each method
}
\description{
\code{spec_molNetwork} uses \code{Spectra} objects and calculates molecular similarity
networks based on the comparsion method(s) defined in \code{methods}
}
\author{
Michael Witting
}
