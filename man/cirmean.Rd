\name{cirmean}
\alias{cirmean}

\title{Circular Mean}

\description{This function computes the circular mean as defined in Mardia et al. (2000).
\deqn{\overline{\theta}=\left\{\begin{array}{ll} \arctan{\Big(\frac{\overline{S}}{\overline{C}}}\Big) \hspace{1cm} &\textrm{si  } \overline{S}>0, \overline{C}>0\\
\arctan{\Big(\frac{\overline{S}}{\overline{C}}}\Big)+\pi \hspace{1cm} &\textrm{si  } \overline{C}<0\\
\arctan{\Big(\frac{\overline{S}}{\overline{C}}}\Big)+2\pi \hspace{1cm} &\textrm{si  } \overline{S}<0, \overline{C}>0\\
\end{array}\right.}{}
\deqn{where\hspace{0.5cm} \overline{S}=\frac{1}{n}\sum_{i=1}^{n}\sin{\theta_{i}}\hspace{0.5cm} and \hspace{0.5cm} \overline{C}=\frac{1}{n}\sum_{i=1}^{n}\cos{\theta_{i}} }{}
}

\usage{cirmean(data)}

\arguments{
  \item{data }{a numeric vector with the circular data.}
}
\details{It is necessary that all the values belong to the interval \eqn{[0,2\pi]}{[0,2pi]}, since the data must belong to the circle.  In this case the arithmetic mean is not a good measure of location as it can be seen in the references.  So the circular mean that this function computes should be used to calculate the mean of angles.
}

\value{The value of the circular mean.}

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
\code{\link{CIREi}}, \code{\link{cirSCE}}, \code{\link{CTi}}, \code{\link{cirkappa}}, \code{\link{mrl}}.
}

\examples{
data(cirdata)
# We calculate the circular mean of these values:
cirmean(cirdata)
# We calculate the circular means of each population when data have replications.
data(datareplic)
apply(datareplic,1,cirmean)
}

\keyword{circular}
\keyword{mean}