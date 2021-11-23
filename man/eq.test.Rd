\name{eq.test}
\alias{eq.test}

\title{Test of Equality of Circular Orders}

\description{
This function calculates the test of equality of circular orders in different populations.
}

\usage{eq.test(data, popu, ws=NULL, method=NULL, control.method=NULL, output=NULL, coef=1, N=500)}

\arguments{
  \item{data}{matrix or data.frame with the data. See details.}
  \item{popu}{a numeric vector with population to each experiment belongs.}
  \item{ws}{a numeric vector with the values to be used as weight per experiment.}
  \item{method}{The method to be used to aggregate circular orders with ACO function.}
\item{control.method}{The argument to control the method in ACO function.}
\item{output}{The path where write the output of the global orders.}
\item{coef}{The coefficient to use in case of method=TSP, by default 1.}
\item{N}{The number of randomization selections, by default 500.}
}

\details{
This function performs the test to constrast equality of circular orders:

\eqn{\hspace{1.7cm}H_{0}:}{Ho:} The circular parameters follow the same circular order in all populations.

\eqn{\hspace{1.7cm}H_{1}:\hspace{0.1cm}H_{0}}{H1: Ho:} is not true.

The \code{data} must have the elements in the columns and the experiments in the rows.
}

\value{
The output is a list with the following values:
\item{allorders}{matrix, in each row the circular order obtained with all selected experiments and the value for the statistic test in that selection.}
\item{pvalue}{numeric, it is the p-value what results of the test.}
\item{global_order}{numeric vector with the elements ordered as the global circular order estimate.}
\item{CC}{numeric, the confidence coefficient (in percentage) of the global order.}
\item{MFO}{numeric vector with the elements ordered as the Most Frequent global Order in the randomization procedure.}
\item{CCMFO}{numeric, the confidence coefficient (in percentage) of the Most Frequent global Order.}
Two additional outputs could be obtained in the form of .csv files written in the path given by the user in the argument \code{output}:
\item{globalorders.csv}{all the global orders obtained from the randomization procedure and the value of the statistic in each selection.}
\item{frequencydist.csv}{the frequency distribution of all the global orders.} 
}


\author{Author(s): Sandra Barrag?n.
Maintainer:<sandra.barragan@gmail.com>}


\references{

Rueda, C., Fernandez, M. A. and Peddada, S. D. (2009).
  Estimation of parameters subject to order restrictions on a circle with application to estimation of phase angles of cell-cycle genes.
  \emph{Journal of the American Statistical Association}, \bold{104}, n485; pp 338--347. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2742472/}

Fernandez, M. A., Rueda, C. and Peddada, S. D. (2012).
    Identification of a core set of signature cell cycle genes whose relative order of time to peak expression is conserved across
species, \emph{Nucl. Acids Res.} \bold{40}, n7: pp 2823--2832. doi:10.1093/nar/gkr1077.
\url{https://academic.oup.com/nar/article/40/7/2823/1183140}

Barragan, S., Rueda, C., Fernandez, M.A. and Peddada, S.D. (2015).
    Determination of Temporal Order among the Components of an Oscillatory System.
    \emph{PLOS ONE}. \bold{10}, n7: pp 1--14. doi: 10.1371/journal.pone.0124842. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4495067/}

}

\seealso{
\code{\link{ACO}}, \code{\link{CIRE}},\code{\link{sce}}, \code{\link{mrl}}, \code{\link{isocir}}, \code{\link{plot.isocir}}.
}


\examples{
data(cirgenes)
eq.test(cirgenes[-8,c(1:5)], popu=c(rep(1,5),rep(2,4)),
ws=c(1,2,3.5,2,1,8,4.2,1.35,6), method="TSP",control.method="alpha3",N=2)

}

\keyword{circular}
\keyword{isotonic}
\keyword{CIRE}
\keyword{test}
\keyword{equality}
\keyword{order}
\keyword{aggregation}