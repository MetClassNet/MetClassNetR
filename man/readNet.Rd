% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NetworkComparisonFunctions.R
\name{readNet}
\alias{readNet}
\title{Read all the networks stored in a given directory}
\usage{
readNet(netDir, directed = FALSE, pattern = "", format = "csv")
}
\arguments{
\item{netDir}{\code{Path}, directory containing all the networks. NOTE. The networks must be
stored in csv format}

\item{directed}{\code{boolean}, value indicating if the networks are directed or not, the default
value is FALSE (i.e., undirected networks)}

\item{pattern}{\code{string}, (optional) pattern that the file names must contain to be read.
NOTE. If all the networks in \code{netDir} are to be read, the parameter
\code{pattern} can be omitted}

\item{format}{\code{string}, files' format (i.e., extension of the network files). NOTE. If the
format is "csv", then the files must contain 2 columns (source - target) and
each row must be an edge}
}
\value{
List of igraph objects, one per network to compare
}
\description{
The function \code{readNet} reads all the networks stored in a given directory
}
\examples{
# See the NetworkComparison vignette

}
\author{
Elva Novoa, \email{elva-maria.novoa-del-toro@inrae.fr}
}
