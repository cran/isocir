\name{msce}
\alias{msce}
\title{Mean Sum of Circular Errors.}

\description{
This function computes the MSCE between and order and a set of data.
}

\usage{
msce(data, posorder, ws=NULL, ...)
}

\arguments{
  \item{data}{vector or matrix of data to be processed.  See details.}
  \item{posorder}{The posicions of the order.}
\item{...}{arguments to pass to CIRE.}
\item{ws}{the weights.}
}

\details{
This function computes the MSCE between and order and a set of data.

}

\value{
A list with the elements: 
 \item{msce}{The mean of the mean sum of circular errors.}
 \item{msces}{a vector with the mean sum of circular errors for each experiment.}
}


\author{Author(s):  Sandra Barragán.  Maintainer:  <sandraba@eio.uva.es>
}


\seealso{
 \code{\link{sce}}, \code{\link{cond.test}}, \code{\link{mrl}}, \code{\link{isocir}}, \code{\link{plot.isocir}}.
}

\examples{
data(cirgenes)
msce(cirgenes, posorder = c(1:ncol(cirgenes)))
}

\keyword{circular}
\keyword{isotonic}
\keyword{order}
\keyword{CIRE}