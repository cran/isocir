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
\deqn{\overline{R}=\frac{1}{n}\sqrt{S^{2}+C^{2}}}{R_bar=(1/n)*sqrt(S_^2+C^2)}
\deqn{where\hspace{0.5cm} S=\sum_{k=1}^{n}\sin{\theta_k}\hspace{0.5cm} and \hspace{0.5cm} C=\sum_{k=1}^{n}\cos{\theta_k}}{S=sum_{k=1}^{n}(sin(theta_k)) and C=sum_{k=1}^{n}(cos(theta_k))}

The argument \code{data} could be a matrix with \emph{n} columns and \emph{q} rows, \emph{q} is the number of populations.  \code{data} could also be a vector.  For both cases the function \code{\link{rho.circular}} from the package \bold{circular} is used in the calculations.

Missing values in the replications are allowed.
}
\value{
\item{mrl}{numeric vector of dimension \emph{q} with the mean resultant lengths. The \emph{i} element is the mean resultant lenght of the \emph{i} population which is in the row \emph{i} of the matrix \code{data}.}
}

\author{Author(s):  Sandra Barrag?n.  Maintainer:  <sandraba@eio.uva.es>
}

\references{

Mardia, K. and Jupp, P. (2000).
   \emph{Directional Statistics}, Chichester: Wiley.


Rueda, C., Fernandez, M. A. and Peddada, S. D. (2009).
  Estimation of parameters subject to order restrictions on a circle with application to estimation of phase angles of cell-cycle genes.
  \emph{Journal of the American Statistical Association}, \bold{104}, n485; pp 338--347. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2742472/}

Fernandez, M. A., Rueda, C. and Peddada, S. D. (2012).
    Identification of a core set of signature cell cycle genes whose relative order of time to peak expression is conserved across
species, \emph{Nucl. Acids Res.} \bold{40}, n7: pp 2823--2832. doi:10.1093/nar/gkr1077.
\url{http://nar.oxfordjournals.org/content/early/2011/11/30/nar.gkr1077}

}

\seealso{
\code{\link{CIRE}}, \code{\link{sce}}, \code{\link{cond.test}}.
}

\examples{
data(datareplic)
mrl(datareplic)
}

\keyword{circular}
\keyword{resultant}
\keyword{length}
\keyword{CIRE}
