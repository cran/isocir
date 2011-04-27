\name{mrl}
\alias{mrl}

\title{Mean Resultant Length}

\description{
This function calculates the mean resultant length as defined in Mardia et al. (2000). 
}

\usage{mrl(data)}

\arguments{
  \item{data}{matrix with the data}
 }

\details{
It is supposed that we have \emph{n} replications for each population.
\deqn{\overline{R}=\sqrt{S^{2}+C^{2}}}{R=sqrt(S_^2+C^2)}
\deqn{where\hspace{0.5cm} S=\sum_{k=1}^{n}\sin{\theta_k}\hspace{0.5cm} and \hspace{0.5cm} C=\sum_{k=1}^{n}\cos{\theta_k}}{S=sum_{k=1}^{n}(sin(theta_k)) and C=sum_{k=1}^{n}(cos(theta_k))}

The argument \code{data} is a matrix with \emph{n} columns and \emph{q} rows, \emph{q} is the number of populations.

Missing values in the replications are allowed.
}
\value{
\item{mrl}{numeric vector of dimension \emph{q} with the mean resultant lengths. The \emph{i} element is the mean resultant lenght of the \emph{i} population which is in the row \emph{i} of the matrix \code{data}.}
}

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
\code{\link{CIREi}}, \code{\link{cirmean}}, \code{\link{cirSCE}}, \code{\link{CTi}}, \code{\link{cirkappa}}.
}

\examples{
data(datareplic)
mrl(datareplic)
}

\keyword{circular}
\keyword{resultant}
\keyword{length}
\keyword{CIRE}
