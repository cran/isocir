\name{cond.test}
\alias{cond.test}

\title{Conditional Test for Contrasting Circular Order}

\description{
This function calculates the p-value corresponding to the conditional test where the hypotheses are whether the data follow a fixed circular order or not.
}

\usage{cond.test(data, groups=c(1:nrow(data)), kappa=NULL, bias=TRUE)}

\arguments{
  \item{data}{vector or matrix with the data. See details.}
  \item{groups}{a numeric vector with the order to be contrasted.}
  \item{kappa}{numeric and positive, if NULL, it will be estimated.}
  \item{bias}{logical, if \code{TRUE}, the correction of the bias is done in the estimation of kappa.}
}

\details{
This function performs the conditional test to constrast the hypotheses:

\eqn{\hspace{1.7cm}H_{0}:}{Ho:} The circular parameters follow a fixed circular order.

\eqn{\hspace{1.7cm}H_{1}:\hspace{0.1cm}H_{0}}{H1: Ho:} is not true.

This conditional test is a modification of the likelihood ratio test.  See Fernandez et al. (2011) for full details about the conditional test.  For a whole explanation of the use of this function, see the vignette of this package where some examples are shown.

The result of \code{cond.test} is the p-value of this test.  There are two cases depending on whether \eqn{\kappa}{kappa} is known or not.  

If \eqn{\kappa}{kappa} is known, \code{data} is usually a vector with the unrestricted circular means. In this case the value of \eqn{\kappa}{kappa} has to be introduced in the argument \code{kappa}.  \code{data} could also be a matrix with replications and the user knows the value of \eqn{\kappa}{kappa}, then if it is introduced in the argument \code{kappa} it will be used by the function to perform the contrast instead of using the estimation.

When \eqn{\kappa}{kappa} is unknown, replications are needed in order to calculate the estimation. Then, \code{data} must be a matrix where each column is a replication while each row is a population.  \eqn{\kappa}{kappa} is estimated internally using the function \code{mle.vonmises} developed in the package \bold{circular}.  By default, the correction of the bias is done, if the user prefers not correcting the bias, just set the argument \code{bias=FALSE}.

For both cases (\eqn{\kappa}{kappa} known or unknown), the order to be contrasted is introduced in \code{groups} in the form of a numeric vector with length equal to the number of populations. Each value of the argument establishes the level of the order the corresponding population is assigned to.  The default value is the sequence corresponding to the simple order (1, 2, ..., q).  Notice that no particular order is assumed among the populations belonging to the same level set. 
}

\value{
The output is an S3 object of class \code{isocir}, similar to the result in the function \code{\link{CIRE}} but adding the following values:
\item{pvalue}{numeric, it is the p-value what results of the conditional test.}
\item{kappa}{numeric, it is the value of kappa.}
The attribute \code{estkappa} shows if kappa is estimated or known.
}


\author{Author(s): Sandra Barrag�n based on the SAS code written by Miguel A. Fern�ndez.
Maintainer:<sandraba@eio.uva.es>}


\references{

Mardia, K. and Jupp, P. (2000).
   \emph{Directional Statistics}, Chichester: Wiley.

Rueda, C., Fern�ndez, M. A. and Peddada, S. D. (2009).
  Estimation of parameters subject to order restrictions on a circle with application to estimation of phase angles of cell-cycle genes.
  \emph{Journal of the American Statistical Association}, \bold{104}, n485; pp 338--347.
  \url{http://pubs.amstat.org/doi/abs/10.1198/jasa.2009.0120}

Fern�ndez, M. A., Rueda, C. and Peddada, S. D. (2011).
    Identification of a core set of signature cell cycle genes whose relative order of time to peak expression is conserved across
species, \emph{Nucl. Acids Res.} \bold{40}, n7: pp 2823--2832. doi:10.1093/nar/gkr1077.
\url{http://nar.oxfordjournals.org/content/early/2011/11/30/nar.gkr1077}
}

\seealso{
\code{\link{CIRE}},\code{\link{sce}}, \code{\link{mrl}}, \code{\link{isocir}}, \code{\link{plot.isocir}}.
}


\examples{
data(cirdata)
# Example without replications and a general isotropic order:
cond.test(cirdata, groups=c(1,2,1,3,3,4,5,6), kappa=0.2)
# Example with replications and the isotropic order (by default):
data(datareplic)
cond.test(data=datareplic)
}

\keyword{circular}
\keyword{isotonic}
\keyword{CIRE}
\keyword{test}