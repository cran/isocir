\name{CIREi}
\alias{CIREi}
\title{Circular Isotonic Regresssion Estimator}

\description{
This function calculates the Circular Isotonic Regression Estimator (CIRE) of a set of data assuming an order among them.  See details for a full explanation.
}

\usage{
CIREi (data, levels=c(1:nrow(data)), isotropic = TRUE, 
       graphic=FALSE, stack=TRUE, ...)
}

\arguments{
  \item{data}{a matrix of data to be processed.  See details.}
  \item{levels}{a numeric vector with the assumed order of data.}
  \item{isotropic }{logical,  by default \code{TRUE} is the order closed in the circle (isotropic order). If \code{FALSE} the order is not closed.}
  \item{graphic}{logical, the plots are drawn or not, by default FALSE.}
  \item{stack}{logical, it stacks the data in the plot with the default value TRUE.}
  \item{\dots}{further parameters for the plots.}
}

\details{
 It uses the algorithm developed in Rueda et al (2009).  The final solution is found in a computationally efficient way by searching for the feasible solution with the smallest SCE, see \code{\link{cirSCE}}.  The order can be chosen by the user and is set using the \code{levels} argument.  A general case would be the general isotropic order, defined by the set \eqn{C_{GIO}}{C_GIO}.

\deqn{ C_{GIO}=\left\{ \phi\in [0,2\pi]^q:\left \{
  \begin{array}{c}
     \phi_{11} \\
    \phi_{12} \\
    \vdots \\
    \phi_{1 l_{1}}  \\
  \end{array}
\right\} \preceq \left \{
  \begin{array}{c}
     \phi_{21} \\
    \phi_{22} \\
    \vdots \\
    \phi_{2 l_{2}}  \\
  \end{array}
\right\} \preceq \ldots \preceq \left\{
  \begin{array}{c}
     \phi_{L1} \\
    \phi_{L2} \\
    \vdots \\
    \phi_{L l_{L}}  \\
  \end{array}
\right\} \preceq \left \{
  \begin{array}{c}
     \phi_{11} \\
    \phi_{12} \\
    \vdots \\
    \phi_{1 l_{1}}  \\
  \end{array}
\right\} \right\} 
 }{C_GIO=\{phi \in [0,2\pi]^q: [phi_{11}, phi_{12},...phi_{1 l_{1}}]<= [phi_{21},phi_{22},...,phi_{2 l_{2}}]<= ...<= [phi_{L1},phi_{L2},....,phi_{L l_{L}}]<= [phi_{11}, phi_{12},...phi_{1 l_{1}}]\}}

where \emph{L} is the number of levels in the order, \eqn{l_j} is the number of angles in the level \emph{j} and \emph{q} is the total number of angles \eqn{: q=\sum_{j}^{L}l_{j}.}{.}

An example of a general non isotropic order would be:
\deqn{ C_{GNIO}=\left\{ \phi\in [0,2\pi]^q:\left \{
  \begin{array}{c}
     \phi_{11} \\
    \phi_{12} \\
    \vdots \\
    \phi_{1 l_{1}}  \\
  \end{array}
\right \} \preceq \left\{
  \begin{array}{c}
     \phi_{21} \\
    \phi_{22} \\
    \vdots \\
    \phi_{2 l_{2}}  \\
  \end{array}
\right\} \preceq \ldots \preceq \left \{
  \begin{array}{c}
     \phi_{L1} \\
    \phi_{L2} \\
    \vdots \\
    \phi_{L l_{L}}  \\
  \end{array}
\right\} \right\} 
 }{C_NIO=\{phi \in [0,2\pi]^q: [phi_{11}, phi_{12},...phi_{1 l_{1}}]<= [phi_{21},phi_{22},...,phi_{2 l_{2}}]<= ...<= [phi_{L1},phi_{L2},....,phi_{L l_{L}}]\}}

For further information see the vignette of this package \bold{isocir}.


The argument \code{data} is a matrix whose columns are the replications of the data, if they exist, and the rows are the populations.  All the values must belong to the interval \eqn{[0,2\pi]}{[0,2pi]}.  

The argument \code{levels} is a numeric vector with length equal to the number of populations, \emph{q}. Each value of the argument establishes the level of the order the corresponding population is assigned to.  The default value is the sequence corresponding to the simple order (1, 2, ..., q).  Notice that no particular order is assumed among the populations belonging to the same level set.

}

\value{
  \item{cirmeans}{A list with the circular means in the order introduced by \code{levels}.}
  \item{SCE}{Sum of Circular Error of the final CIRE.}
  \item{CIRE}{A list with the ordered estimations, the result of the algorithm that obtained the CIRE.}
}

\author{Author(s):  The implementation in R is done by Sandra \enc{Barragan}{Barragan} based on the SAS routines written by Miguel A. \enc{Fernandez}{Fernandez}.  Maintainer:  <sandraba@eio.uva.es>
}

\references{
Mardia, K. and Jupp, P. (2000).
   \emph{Directional Statistics}, Chichester: Wiley.

Rueda, C., Fernandez, M. A. and Peddada, S. D. (2009).
  Estimation of parameters subject to order restrictions on a circle with application to estimation of phase angles of cell-cycle genes.
  \emph{Journal of the American Statistical Association}, \bold{104}, n485; pp 338--347.
  \url{http://pubs.amstat.org/doi/abs/10.1198/jasa.2009.0120}

\enc{Fernandez}{Fernandez}, M. A., Rueda, C. and Peddada, S. D. (2011).
   Isotropic order among core set of orthologs conserved between budding and fission yeasts. \emph{Prepint}.
}

\seealso{
\code{\link{cirmean}}, \code{\link{cirSCE}}, \code{\link{CTi}}, \code{\link{cirkappa}}, \code{\link{mrl}}.
}

\examples{

# We consider the following data from the package that are random circular data:
data(cirdata)
data1<-cbind(cirdata)
circular_ordered_estimator <- CIREi(data=data1)
# We can take the vector of the CIRE estimators:
circular_ordered_estimator $CIRE
# And the SCE:
circular_ordered_estimator $SCE

# Random data with a more complex order:
CIREi(data1, levels=c(1,1,2,3,5,3,4,6))

}

\keyword{circular}
\keyword{isotonic}
\keyword{isotropic}
\keyword{CIRE}
