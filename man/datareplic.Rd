\name{datareplic}
\alias{datareplic}

\docType{data}
\title{Random Circular Data with Replications.}

\description{An example data set.   The data have been created randomly for the purpose of running some examples of the functions of the package \bold{isocir}.}

\usage{data(datareplic)}

\details{
  Every data of the matrix come from the unit circle \eqn{[0,2\pi]}{[0,2pi]}.  

The data in the matrix \code{datareplic} are random angular measurements.  Its columns are the 10 replications and its rows are the 8 populations.
}

\references{
 
Mardia, K. and Jupp, P. (2000).
  \emph{Directional Statistics}, Chichester: Wiley.

}

\seealso{
\code{\link{cirgenes}}, \code{\link{cirdata}}.
}


\examples{
data(datareplic)

# example testing the next order:
orderGroups2 <- c(1:8)

example2test <- cond.test(datareplic, groups = orderGroups2, bias = TRUE)
example2test


round(unlist(example2test$cirmeans), digits = 3)

}

\keyword{datasets}
\keyword{circular}
