\name{CLM}
\alias{CLM}
\title{Circular Local Minimization}

\description{
This function computes the Circular Local Minimization.
}

\usage{
CLM(data, order0, ws=NULL, control.method=c("msce","cirktau"))
}

\arguments{
  \item{data}{matrix of circular data to be processed.}
  \item{order0}{the initial order to be improved.}
\item{ws}{the weights.}
\item{control.method}{by default "msce", maximization of "cirktau".}
}

\details{
This function computes the Circular Local Minimization of the MSCE by default. It computes the Circular Local Maximization of the Circular Kendall Tau if \code{control.method="cirktau"}.

}

\value{
This function returns a list with the next elements:
\item{order0}{the initial order introduced in the arguments.}
\item{msce0}{the mean sum of circular errors of the initial order with the data.}
or
item{tau0}{the mean circular Kendall Tau of the initial order with the data.}
\item{final_order}{the final order after proccess the Circular Local Minimization.}
\item{bestsce}{the msce of the final order with the data.}
or
\item{bestau}{the mean circular Kendall Tau of the final order with the data.}
}


\author{Author(s):  Sandra Barrag?n.  Maintainer:  <sandraba@eio.uva.es>
}

\references{
DWORK, C., KUMAR, R., NAOR, M. and SIVAKUMAR, D. (2001), Rank Aggregation Methods for the
Web, Proceedings of the 10th International World Wide Web Conference, pp. 613--622.
}

\seealso{
 \code{\link{sce}}, \code{\link{cond.test}}, \code{\link{mrl}}, \code{\link{isocir}}, \code{\link{plot.isocir}}.
}

\examples{
data(cirgenes)
aggre_order <- ACO(cirgenes[,c(1:5)], method="TSP", control.method="time")$TSP_order[1,]
CLM(cirgenes[,c(1:5)], order0 = c(1:5))
#datos<- rbind (c(0, 1/10, 1/9, 1, 11/10, 10/9)*pi,c(0, 1/2, 1/10, 1, 3/2, 11/10)*pi) 
#CLM(datos, order0 = c(1:6))
}

\keyword{circular}
\keyword{isotonic}
\keyword{order}
\keyword{CLM}
\keyword{aggregation}
