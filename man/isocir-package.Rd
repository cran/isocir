\name{isocir-package}
\alias{isocir-package}
\docType{package}
\title{
 Isotonic Inference with Circular data.
}
\description{
  
 This is an R package to make ISOtonic inference for CIRcular data.  
}
\details{

\tabular{ll}{
Package: \tab isocir\cr
Version: \tab 1.0\cr
Date: \tab 2011-04-15\cr
Depends: \tab combinat, circular\cr
License: \tab GPL (>=3)\cr
}
%  The functions are:
%  \tabular{ll}{
%  \code{\link{CIREi}} \tab computes the algorithm to obtain the CIRE.\cr
%  \code{\link{CTi}} \tab computes the conditional test for isotropic order.\cr
%  \code{\link{cirmean}} \tab computes the Circular Mean.\cr
%  \code{\link{cirSCE}} \tab computes the Sum of Circular Error.\cr
%  \code{\link{cirkappa}} \tab computes the estimation of kappa.\cr
%  \code{\link{mrl}} \tab computes the estimation of the Mean Resultant Length.\cr
%  }
%  The data sets are:
%  \tabular{ll}{
%  \code{\link{rancir}} \tab set of Random Circular Data.\cr
%  \code{\link{cirgenes}} \tab a matrix, called experiments.\cr
%   }
}
\author{Author(s):  The implementation in R is done by Sandra Barragan based on the SAS code written by Miguel A. Fernandez.  Maintainer:  <sandraba@eio.uva.es>
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
\code{\link{CIREi}}, \code{\link{cirmean}}, \code{\link{cirSCE}}, \code{\link{CTi}}, \code{\link{cirkappa}}, \code{\link{mrl}}, \code{\link{datareplic}}, \code{\link{cirdata}}, \code{\link{cirgenes}}.
}

\keyword{ package }
\keyword{ isocir }