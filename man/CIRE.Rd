\name{CIRE}
\alias{CIRE}
\title{Circular Isotonic Regresssion Estimator}

\description{
This function obtains the Circular Isotonic Regression Estimator (CIRE) of a set of data assuming an order among them.  See details for a full explanation.
}

\usage{
CIRE (data, groups=c(1:nrow(data)), circular = TRUE)
}

\arguments{
  \item{data}{vector or matrix of data to be processed.  See details.}
  \item{groups}{a numeric vector with the assumed order of data.}
  \item{circular}{logical,  by default \code{TRUE} is the order wrapped around the circle (circular order). If \code{FALSE} the order is not closed.}
}

\details{
 It uses the algorithm developed in Rueda et al (2009).  The final solution is found in a computationally efficient way by searching for the feasible solution with the smallest SCE, see \code{\link{sce}}.  The order can be chosen by the user and is set using the \code{groups} argument.  An example of order would be the simple circular order, defined by the set \eqn{C_{sco}}{C_sco}.

\deqn{ C_{sco}=\{\phi=(\phi_1,\phi_2,\ldots,\phi_q)\in
[0,2\pi]^q:\phi_{1}\leq\phi_{2}\leq\ldots\leq\phi_{q}\leq\phi_{1}\}
 }{C_sco=\{phi \in [0,2\pi]^q: phi_1 <= phi_2 <= ...<= phi_q <= phi_1\}}

This order is wrapped around the circle.  However, if the user wants the order not to be wrapped, then \code{circular = FALSE} and the order would be the simple order, for instance:
\deqn{ C_{so}=\{\phi=(\phi_1,\phi_2,\ldots,\phi_q)\in
[0,2\pi]^q:0\leq\phi_{1}\leq\phi_{2}\leq\ldots\leq\phi_{q}\leq 2\pi\}
 }{C_so=\{phi \in [0,2\pi]^q: 0 <= phi_1 <= phi_2 <= ...<= phi_q <= 2pi\}}

For further information see the vignette of this package \bold{isocir}.

The argument \code{data} is a vector with the unrestricted circular means or a matrix whose columns are the replications of the data and the rows are the populations.  All the values must belong to the interval \eqn{[0,2\pi]}{[0,2pi]}.  

The argument \code{groups} is a numeric vector with length equal to the number of populations, \emph{q}. Each value of the argument establishes the level of the order the corresponding population is assigned to.  The default value is the sequence corresponding to the simple order (1, 2, ..., q), circular if \code{closed=TRUE}.  Notice that no particular order is assumed among the populations belonging to the same level set.

}

\value{
It returns an S3 object of class \code{isocir} that is a list with the next elements:
  \item{cirmeans}{A list with the unrestricted circular means in the same form as the order introduced by \code{groups}.}
  \item{SCE}{Sum of Circular Errors of the final CIRE.}
  \item{CIRE}{A list with the ordered estimations, the result of the algorithm that obtained the CIRE.}
The Circular Isotonic Regression Estimator is returned as a list just to see the groups of the order. In case you would like to use it as a vector just run: \code{unlist(object$CIRE)}.
}


\author{Author(s):  Sandra Barragán based on the SAS routines written by Miguel A. Fernández.  Maintainer:  <sandraba@eio.uva.es>
}

\references{
Mardia, K. and Jupp, P. (2000).
   \emph{Directional Statistics}, Chichester: Wiley.

Rueda, C., Fernández, M. A. and Peddada, S. D. (2009).
  Estimation of parameters subject to order restrictions on a circle with application to estimation of phase angles of cell-cycle genes.
  \emph{Journal of the American Statistical Association}, \bold{104}, n485; pp 338--347.
  \url{http://pubs.amstat.org/doi/abs/10.1198/jasa.2009.0120}

Fernández, M. A., Rueda, C. and Peddada, S. D. (2011).
   Identification of a core set of signature cell cycle genes whose relative order of time to peak expression is conserved across
species, \emph{Nucl. Acids Res.} \bold{40}, n7: pp 2823--2832. doi:10.1093/nar/gkr1077.
\url{http://nar.oxfordjournals.org/content/early/2011/11/30/nar.gkr1077}
}

\seealso{
 \code{\link{sce}}, \code{\link{cond.test}}, \code{\link{mrl}}, \code{\link{isocir}}, \code{\link{plot.isocir}}.
}

\examples{

# We consider the following data from the package that are random circular data:
data(cirdata)
circular_ordered_estimator <- CIRE(cirdata)
# We can take the vector of the CIRE estimators:
circular_ordered_estimator $CIRE
# And the SCE:
circular_ordered_estimator $SCE

# Random data with a more complex order:
CIRE(cirdata, groups=c(1,1,2,3,5,3,4,6))

}

\keyword{circular}
\keyword{isotonic}
\keyword{order}
\keyword{CIRE}
