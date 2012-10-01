\name{cirdata}
\alias{cirdata}

\docType{data}
\title{Random Circular Data.}

\description{An example data set.   The data have been created randomly for the purpose of running some examples of the functions of the package \bold{isocir}.}

\usage{data(cirdata)}

\details{
  Every data of the vector come from the unit circle \eqn{[0,2\pi]}{[0,2pi]}.  
}

\references{
 
Mardia, K. and Jupp, P. (2000).
  \emph{Directional Statistics}, Chichester: Wiley.

}

\seealso{
\code{\link{cirgenes}}, \code{\link{datareplic}}.
}


\examples{
data(cirdata)

# Put groups of the order in the correct form:
orderGroups <- c(1,1,1,2,2,3,4,4)


# execute CIRE to obtain the CIRE with circular order:
example1CIRE <- CIRE(cirdata, groups = orderGroups, circular = TRUE)
example1CIRE 

# if the user wants the graphics: 
#plot(example1CIRE)
#plot(example1CIRE, option="cirmeans")



}

\keyword{datasets}
\keyword{circular}
