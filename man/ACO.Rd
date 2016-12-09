\name{ACO}
\alias{ACO}
\title{Aggregation of Circular Orders}

\description{
This function computes the aggregation of circular orders.
}

\usage{
ACO(data, method=c("Naive", "CB", "CMC", "TSP", "CH"), 
control.method, ws=NULL, coef=1)
}

\arguments{
  \item{data}{vector or matrix of data to be processed.  See details.}
  \item{method}{"Naive", "CB", "CMC", "TSP", "CH".}
\item{control.method}{"tau", "MSCE", "pos", "cirmean", "cirmed", 
"1", "2", "3", "4m", "4c", "bin", "pos", "alpha1", "alpha2", "alpha3", "alpha4","alphainf", "time", "arc", "chord", "bin", "pos", "cos", "cmean", "mrl", "e3", "ave", "qua", "nat", "natp", "natb".}
\item{ws}{the weights.}
\item{coef}{The coeficient for TSP method.}
}

\details{
This function computes the aggregation of circular orders.  

The possible methods are:
\code{Naive}: the naive method of choosing the order from the data.
\code{CB}: the circular Borda method.
\code{TSP}: the Traveling Salesperson Problem apply to aggregation of circular orders.
\code{CH}: the circular Hodge method.

The \code{coef} argument is the proportion of maximum presolutions 
check respect to the number of elements n. It means that when the TSP is solved by heuristics, a maximum of the best coef*n 
possible solutions are proned to be the final solution, they are checked and the one with the best SCE is chosen.
}


\value{
A list with the elements: 
\item{aggre_order}{the circular aggregated order.}
\item{msce}{the MSCE corresponding to the circular aggregated order.}
\item{mtau}{the circular Kendall tau corresponding to the circular aggregated order.}
\item{mintour}{the order corresponding to the tour of minimum length in case of TSP.}
\item{mt_msce}{in case of TSP, the MSCE corresponding to the tour of minimum length.}
\item{tour_length}{in case of TSP, the length of the tour.}
\item{scores}{in case of TSP, the length of each edge. In case of circular Hodge, the scores.}
\item{out}{in case of circular Hodge, the element where the circle is cut.}
\item{error1}{in case of circular Hodge method, the error1.}
\item{error2}{in case of circular Hodge method, the error2.}
}


\author{Author(s):  Sandra Barrag?n.  Maintainer:  <sandraba@eio.uva.es>
}

\references{
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
\code{\link{eq.test}}, \code{\link{sce}}, \code{\link{cond.test}}, \code{\link{mrl}}, \code{\link{isocir}}, \code{\link{plot.isocir}}.
}

\examples{
data(cirgenes)
ACO(cirgenes[,c(1:5)], method="TSP", control.method="alpha3")
#datos<- rbind (c(0, 1/10, 1/9, 1, 11/10, 10/9)*pi,c(0, 1/2, 1/10, 1, 3/2, 11/10)*pi) 
#ACO(datos, method="TSP", control.method="alpha3")

}

\keyword{circular}
\keyword{isotonic}
\keyword{order}
\keyword{CIRE}
\keyword{aggregation}
