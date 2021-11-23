\name{plot.isocir}
\alias{plot.isocir}

\title{S3 Method to Plot S3 Objects of Class isocir}
\description{
This is a function to plot an object of class \code{isocir} that is the result of the function \code{\link{CIRE}} and \code{\link{cond.test}}.  It creates a plot of circular data points on the current graphics device.  Data points are either plotted as points on the unit circle, or the range of the circle is divided into a specified number of bins, and points are stacked in the bins corresponding to the number of observations in each bin. 
}

\usage{
\method{plot}{isocir}(x, option=c("CIRE", "cirmeans"), cex = 1, stack = TRUE,
  axes = TRUE, sep = 0.025, shrink = 1, bins = 300, ticks = FALSE,
  tcl = 0.025, tcl.text = 0.125, col = NULL, tol = 0.04, uin = NULL,
  xlim = c(-1, 1), ylim = c(-1, 1), digits = 2, units = NULL,
  template = NULL, zero = NULL, rotation = NULL, 
  main = NULL, sub=NULL, xlab = "", ylab = "", 
  control.circle=circle.control(), ...)
}

\arguments{
\item{x}{an object of class \code{isocir}.}
\item{option}{The points that are used in the plot: If \code{CIRE}, the Circular Isotonic Regression Estimator. If \code{cirmeans}, the unrestricted circular means. Default is \code{CIRE}.}
\item{cex}{point character size.  See help on \code{\link{par}}.}
\item{stack}{logical; if \code{TRUE}, points are stacked on the perimeter of the circle.  Otherwise, all points are plotted on the perimeter of the circle.  Default is \code{TRUE}.}
\item{axes}{logical; if \code{TRUE} axes are plotted according to properties of \code{x}.}
\item{sep}{constant used to specify the distance between stacked points,
  if \code{stack==TRUE} or in the case of more than one dataset.  Default is 0.025; smaller values will create smaller spaces.}
\item{shrink}{parameter that controls the size of the plotted circle.
  Default is 1.  Larger values shrink the circle, while smaller values
  enlarge the circle.}
\item{bins}{if \code{stack==TRUE}, bins is the number of arcs to partition the circle with. Default is 300.}
\item{ticks}{logical; if \code{TRUE} ticks are plotted according to the
  value of \code{bins}.}
\item{tcl}{length of the ticks.}
\item{tcl.text}{The position of the axis labels.}
\item{col}{color of the points. The values are recycled if needed.}
\item{tol}{proportion of white space at the margins of plot.}
\item{uin}{desired values for the units per inch parameter. If of length
  1, the desired units per inch on the x axis.}
\item{xlim, ylim}{the ranges to be encompassed by the x and y
  axes. Useful for centering the plot.}
\item{digits}{number of digits used to print axis values.}
\item{main, sub, xlab, ylab}{title, subtitle, x label and y label of the plot.}
\item{units}{the units used in the plot.}
\item{template}{the template used in the plot.}
\item{zero}{the zero used in the plot.}
\item{rotation}{the rotation used in the plot.}
\item{control.circle}{parameters passed to \code{\link{plot.default}} in order to draw the circle. The function \code{circle.control} is used to set the parameters.}
\item{\dots}{futher parameters passed to \code{\link{points.default}}.}
}

\details{
When there are many closely distributed observations, stacking is
recommended.  When stacking the points, if there are many points in a particular bin, it may be necessary to shrink the plot of the circle so that all points fit. This is controlled with the parameter \code{shrink}. Generally the parameter \code{sep} does not need adjustment, however, when shrinking the plot, or for a very large number of observations, it may be helpful.
}

\value{
A list with information on the plot: zero, rotation and next.points.
}

\note{Some codes from \code{\link{eqscplot}} in MASS is used besides some from \code{\link{plot.circular}} in the circular package.}

\author{Sandra Barragán based on the code of Claudio Agostinelli and Ulric Lund from the package \bold{circular}.}

\seealso{\code{\link{cond.test}}, \code{\link{CIRE}}, \code{\link{isocir}}.}


\examples{

data(cirdata)
result<-CIRE(cirdata)
plot(result)
plot(result,option="cirmeans")
}

\keyword{circular}
\keyword{plot}
\keyword{CIRE}