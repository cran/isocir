\name{cirKendall}
\alias{cirKendall}
\title{Circular Kendall's Tau}

\description{
This function obtains the circular Kendall's tau between two sets of data.
}

\usage{
cirKendall(phi1, phi2, test = FALSE, control.test = c("noteq", "upper", "lower"))
}

\arguments{
  \item{phi1}{a numeric vector of data to be processed.  See details.}
  \item{phi2}{a numeric vector of data.}
  \item{test}{logical,  by default \code{FALSE}.}
  \item{control.test}{"noteq", "upper" or "lower".}
}

\details{
This function obtains the circular Kendall's tau between two data set.
}

\value{
It returns the value of the circular Kendall's tau.
If \code{test=TRUE}, it returns also the pvalue of the test.
}


\author{Author(s):  Sandra Barragán.  Maintainer:  <sandraba@eio.uva.es>
}

\references{
Fisher,N.I. (1993).\emph{Statistical analysis of circular data.} Cambridge University Press.
}

\seealso{
 \code{\link{sce}}, \code{\link{cond.test}}, \code{\link{mrl}}, \code{\link{isocir}}, \code{\link{plot.isocir}}.
}

\examples{

data(cirgenes)
cirKendall(cirgenes[1,],cirgenes[2,])
# Even with missing values:
cirKendall(cirgenes[1,],cirgenes[8,])
}

\keyword{circular}
\keyword{isotonic}
\keyword{order}

