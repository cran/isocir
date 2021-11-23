\name{mcirktau}
\alias{mcirktau}
\title{Mean Circular Kendall's Tau.}

\description{
This function computes the mean Circular Kendall's tau between a circular order and a data set.
}

\usage{
mcirktau(data, posorder, ws=NULL)
}

\arguments{
  \item{data}{vector or matrix of data to be processed.  See details.}
  \item{posorder}{The posicions of the order.}
\item{ws}{the weights.}
}

\details{
This function computes the mean Circular Kendall's tau between a circular order and a data set.

}

\value{
A list with the elements: 
 \item{mtau}{mean circular Kendall's tau.}
 \item{ntaus}{a vector with the circular Kendall's tau of all the experiment.}
}


\author{Author(s):  Sandra Barragán.  Maintainer:  <sandra.barragan@gmail.com>
}


\seealso{
 \code{\link{sce}}, \code{\link{cond.test}}, \code{\link{cirKendall}}, \code{\link{isocir}}, \code{\link{plot.isocir}}.
}

\examples{
data(cirgenes)
mcirktau(cirgenes, posorder = c(1:ncol(cirgenes)))
}

\keyword{circular}
\keyword{isotonic}
\keyword{order}
\keyword{CIRE}
\keyword{aggregation}