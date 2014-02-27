\name{isocir}
\alias{isocir}
\alias{is.isocir}
\alias{print.isocir}
\title{S3 Objects of Class isocir.}
\description{
This function creates an S3 object of class \code{isocir}.  There are also some methods described for the class \code{isocir}. 
}
\usage{
isocir(cirmeans = NULL, SCE =NULL, CIRE = NULL, pvalue = NULL, kappa = NULL)
is.isocir(x)
\method{print}{isocir}(x, decCIRE = 3, decpvalue = 4, deckappa = 2, ...)
}

\arguments{
  \item{cirmeans}{NULL by default, the unrestricted circular means}
\item{SCE}{NULL by default, the Sum of Circular Errors.} 
\item{CIRE}{NULL by default, the Circular Isotonic Regression Estimator.}
\item{pvalue}{NULL by default, the result of the conditional test}
\item{kappa}{NULL by default.}
  \item{x}{an object of class CIRE.}
  \item{decCIRE}{the number of decimal places to round the CIRE and the SCE.}
  \item{decpvalue}{the number of decimal places to round the pvalue.}
  \item{deckappa}{the number of decimal places to round kappa.}
\item{\dots}{futher parameters passed to \code{print.default}.}

 }

\details{
With the function \code{isocir} a list of class \code{isocir} is created with the elements: \code{cirmeans}, \code{SCE}, \code{CIRE} even if they are NULL.  However, \code{pvalue} and \code{kappa} are added just in the case that they have any value.  

The argument \code{x} must be of class \code{isocir}.
}

\value{
An object of class \code{isocir}. 
}

\author{ Sandra Barragán.}

\seealso{
\code{\link{CIRE}}, \code{\link{cond.test}}.
}

\examples{
data(cirdata)
x <- CIRE(cirdata)
print(x)
is.isocir(x)
class(x)

plot(x)
class(x)

# If you want to use the CIRE in other calculations you can obtain it as a vector:
unlist(x$CIRE)
# But be careful because this unclass and lost attributes! 

# To create a new object of class CIRE:
y <- isocir()
}

\keyword{isocir}
