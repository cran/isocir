\name{cirSCE}
\alias{cirSCE}

\title{Sum of Circular Error}
\description{This function calculates the distance between two \emph{q} dimensional points, or between a point and a set of points, in the circular space.  
}

\usage{cirSCE(arg1, arg2, mrl=1)}
\arguments{
  \item{arg1}{vector with the values for the first point.}
  \item{arg2}{vector (or matrix) with the values for the second point (set of points).}
  \item{mrl}{vector with the mean resultant lenghts, by default 1.}
}

\details{
This function is useful to measure the circular error when an estimator is used to estimate an unknown \emph{q}-dimensional parameter. It plays a role in the circular space similar to that of the mean squared error in the usual Euclidean space.

When \code{arg2} is a vector:  \deqn{SCE(arg1,arg2)=\sum_{i=1}^{q} r_{i} [1-\cos (arg1_{i}-arg2_{i})]}

If \code{arg1} and \code{arg2} are just to \emph{q} dimensional points, both arguments must be vectors with the same length and \code{mrl} is not needed.  

If \code{arg2} is a mean vector coming from several replications, the corresponding mean resultant lengths (\eqn{r_i}) have to be introduced in the argument \code{mrl}. 

When \code{arg2} is a matrix:  \deqn{SCE(arg1,arg2)=\sum_{i=1}^{q}\sum_{k=1}^{n} [1-\cos (arg1_{i}-arg2_{ik})]}

If \code{arg2} is a matrix with the replications in its columns, \code{mrl} is not needed as the mean resultant lengths are calculated internally.   

}

\value{The value of the Sum of Circular Error.}

\author{Author(s):  The implementation in R is done by Sandra \enc{Barragan}{Barragan} based on the SAS routine written by Miguel A. \enc{Fernandez}{Fernandez}.  Maintainer:  <sandraba@eio.uva.es>
}
\references{
Mardia, K. and Jupp, P. (2000).
   \emph{Directional Statistics}, Chichester: Wiley.

Rueda, C., \enc{Fernandez}{Fernandez}, M. A. and Peddada, S. D. (2009).
  Estimation of parameters subject to order restrictions on a circle with application to estimation of phase angles of cell-cycle genes.
  \emph{Journal of the American Statistical Association}, \bold{104}, n485; pp 338--347.
  \url{http://pubs.amstat.org/doi/abs/10.1198/jasa.2009.0120}

\enc{Fernandez}{Fernandez}, M. A., Rueda, C. and Peddada, S. D. (2011).
   Isotropic order among core set of orthologs conserved between budding and fission yeasts. \emph{Prepint}.
}

\seealso{
\code{\link{CIREi}}, \code{\link{cirmean}}, \code{\link{CTi}}, \code{\link{cirkappa}}, \code{\link{mrl}}.
}
\examples{
data(cirdata)
exampledata1 <- cirdata
exampledata2 <- (cirdata+(pi/4))%%(2*pi) #same length
cirSCE(exampledata1,exampledata2)
}

\keyword{SCE}
\keyword{circular}
\keyword{error}
