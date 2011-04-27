\name{cirkappa}

\alias{cirkappa}
\title{Kappa, parameter of the von Mises distribution}

\description{
This function obtains the estimation of kappa, the concentration parameter of a von Mises distribution.
}

\usage{cirkappa(data, biascorr=TRUE)}

\arguments{
  \item{data}{a matrix of data.}
  \item{biascorr}{logical, by default \code{TRUE} when a bias corrected method is applied.}
}

\details{
The estimation is made following Mardia et al. (2000).  It is necessary to have replications in the data set.
\deqn{\kappa=A^{-1}(\overline{R})}{kappa=A^{-1}(R)}
\deqn{where \hspace{0.5cm} \overline{R}=\sqrt{C^{2}+S^{2}},\hspace{0.5cm} with\hspace{0.5cm} S=\sum_{i=1}^{q}\sum_{k=1}^{n}\sin{\theta_{ik}} \hspace{0.5cm} and \hspace{0.5cm} C=\sum_{i=1}^{q}\sum_{k=1}^{n}\cos{\theta_{ik}},}{R=sqrt(C^2+S^2),    with   S=sum_{i=1}^{q}sum_{k=1}^{n}(sin(theta_ik))   and   C=sum_{i=1}^{q}sum_{k=1}^{n}(cos(theta_ik)),}
\deqn{and\hspace{0.5cm}  A(\kappa)=\frac{I_1(\kappa)}{I_0(\kappa)}}{and    A(kappa)=I_1(kappa)/I_0(kappa),}
where \eqn{I_0}{I_0} and \eqn{I_1}{I_1} are the modified Bessel function of the first kind and orders 0 and 1 respectively and \eqn{\theta_{ik}}{\emph{theta_ik}} is the angular observation of replication \emph{k} in population \emph{i}.

The argument \code{data} is a matrix with \emph{n} columns and \emph{q} rows where \emph{n} is the number of replications and \emph{q} is the number of populations.  Obviously, \eqn{N=q*n}{N=q*n}, although missing values in the replications are allowed.  The bias is corrected in case of \code{biascorr=TRUE} following the proposed estimator by Best and Fisher (1986) and appearing in Mardia et al. (2000) [p. 87].

}


\value{
\item{kappa}{numeric, the estimation of the concentration parameter of a von Mises distribution.}

}

\author{Author(s):  The implementation in R is done by Sandra \enc{Barragan}{Barragan} based on the SAS routines written by Miguel A. \enc{Fernandez}{Fernandez}.  Maintainer:  <sandraba@eio.uva.es>
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

Best, D. J. and Fisher, N. I. (1986).
   The bias of the maximum likelihood estimators of the von Mises-Fisher concentration parameters. \emph{Comm. Statist. Simulation Comput.}
   \bold{10}, pp 493--502 (87, 200). 
}

\seealso{
\code{\link{CIREi}}, \code{\link{cirmean}}, \code{\link{cirSCE}}, \code{\link{CTi}}, \code{\link{mrl}}.
}

\examples{
data(datareplic)
cirkappa(datareplic)
}

\keyword{circular}
\keyword{isotonic}
\keyword{vonMises}
\keyword{kappa}
