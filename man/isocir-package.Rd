\name{isocir-package}
\alias{isocir-package}
\docType{package}
\title{
 Isotonic Inference with Circular data.
}
\description{
  
 This is an R package to draw ISOtonic inference for CIRcular data.  
}
\details{

\tabular{ll}{
Package: \tab isocir\cr
Version: \tab 1.9\cr
Depends: \tab combinat, circular\cr
License: \tab GPL (>=3)\cr
}
%  The functions are:
%  \tabular{ll}{
%  \code{\link{CIRE}} \tab computes the algorithm to obtain the CIRE.\cr
%  \code{\link{cond.test}} \tab computes the conditional test for isotropic order.\cr
%  \code{\link{sce}} \tab computes the Sum of Circular Error.\cr
%  \code{\link{mrl}} \tab computes the estimation of the Mean Resultant Length.\cr
%  }
%  The data sets are:
%  \tabular{ll}{
%  \code{\link{datareplic}} \tab set of Random Circular Data with replications.\cr
%  \code{\link{cirdata}} \tab set of Random Circular Data.\cr
%  \code{\link{cirgenes}} \tab a matrix, called experiments.\cr
%   }
}
\author{Author(s):  Sandra Barragán mainly based on the SAS code written by Miguel A. Fernández.  Maintainer:  <sandraba@eio.uva.es>
}
\references{


Mardia, K. and Jupp, P. (2000).
   \emph{Directional Statistics}, Chichester: Wiley.
   
Barragan, S.,  Fernandez, M.A., Rueda, C. and Peddada, S.D. (2015).
isocir: An R Package for Constrained Inference Using Isotonic Regression for Circular Data, with an Application to Cell Biology. \emph{Journal of Statistical Software}, \bold{54}, n1: pp 1--17. doi: 10.18637/jss.v054.i04. \url{https://www.jstatsoft.org/index.php/jss/article/view/v054i04}


Rueda, C., Fernandez, M. A. and Peddada, S. D. (2009).
  Estimation of parameters subject to order restrictions on a circle with application to estimation of phase angles of cell-cycle genes.
  \emph{Journal of the American Statistical Association}, \bold{104}, n485; pp 338--347. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2742472/}

Fernandez, M. A., Rueda, C. and Peddada, S. D. (2012).
    Identification of a core set of signature cell cycle genes whose relative order of time to peak expression is conserved across
species, \emph{Nucl. Acids Res.} \bold{40}, n7: pp 2823--2832. doi:10.1093/nar/gkr1077.
\url{http://nar.oxfordjournals.org/content/early/2011/11/30/nar.gkr1077}

Barragan, S., Rueda, C., Fernandez, M.A. and Peddada, S.D. (2015).
    Determination of Temporal Order among the Components of an Oscillatory System.
    \emph{PLOS ONE}. \bold{10}, n7: pp 1--14. doi: 10.1371/journal.pone.0124842. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4495067/}

}



\seealso{
\code{\link{CIRE}}, \code{\link{sce}}, \code{\link{cond.test}}, \code{\link{mrl}}, \code{\link{datareplic}}, \code{\link{cirdata}}, \code{\link{cirgenes}}.
}

\keyword{ package }
\keyword{ isocir }