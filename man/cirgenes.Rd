\name{cirgenes}
\alias{cirgenes}

\docType{data}
\title{A set of angular measurements from cell-cycle experiments with genes.}

\description{A data set taken from some cell-cycle experiments with the purpose of running some examples of the functions of the package \bold{isocir}.  See details.}

\usage{data(cirgenes)}

\format{
The set has a single element:
  \describe{
    \item{\code{cirgenes}}{a matrix with the measures of several cell-cycle experiments.}
  }
}

\details{
 This data are located in the unit circle \eqn{[0,2\pi]}{[0,2pi]}.  

That matrix \code{cirgenes}  has in its rows the peak expression points for the genes, see Source.  The columns of the matrix are the experiments.  There are 16 genes (rows) and 10 experiments (columns).  The experiments are not held under the same conditions so that the data cannot be considered as replications.
}
\source{
Initially, the data come from the experiments of Oliva et al. (2005), Peng et al. (2005) and Rustici et al. (2004).  They are the observed points of peak expressions of some genes along the cell cycle of the S. Pombe yeast and have been obtained using the Random Periods Model developed in Liu et al. (2004).   The source of that data is http://www.cyclebase.org/CyclebaseSearch, Gauthier (2007).
}
\references{

Fern?ndez, M. A., Rueda, C. and Peddada, S. D. (2012).
    Identification of a core set of signature cell cycle genes whose relative order of time to peak expression is conserved across
species, \emph{Nucl. Acids Res.} \bold{40}, n7: pp 2823--2832. doi:10.1093/nar/gkr1077.
\url{http://nar.oxfordjournals.org/content/early/2011/11/30/nar.gkr1077}

Gauthier, N. (2007). 
Cyclebase.org - A Comprehensive Multi-Organism Online Database of Cell-Cycle Experiments.
\emph{Nucleic Acids Research}. \bold{36}, pp 854--859. http://www.cyclebase.org/CyclebaseSearch.
 
Liu, D., Umbach, D., Peddada, S., Li, L., Crockett, P., Weinberg, C. (2004). 
A Random Periods Model for Expression of Cell-Cycle Genes.
\emph{The National Academy of Sciences of the USA}. \bold{101}(19), pp 7240--7245.

Oliva, A., Rosebrock, A., Ferrezuelo, F., Pyne, S., Chen, H., Skiena, S., Futcher, B. and Leatherwood, J. (2005). 
  The Cell Cycle-Regulated Genes of Schizosaccharomyces pombe.
  \emph{Plos biology}. \bold{3}.

Peng, X., Karuturi, R., Miller, L., Lin, K., Jia, Y., Kondu, P., Wang, L., Wong, L., Liu, E., Balasubramanian, M. and Liu, J. (2005). 
  Identification of Cell Cycle-Regulated Genes in Fission Yeast. 
  \emph{The American Society for Cell Biology}. \bold{16}, pp 1026--1042.

Rustici, G., Mata, J., Kivinen, K., Lio, P., Penkett, C., Burns, G., Hayles, J., Brazma, A., Nurse, P., and Bahler, J. (2004). 
  Periodic Gene Expression program of the fission yeast cell cycle. 
  \emph{Nature Genetics}. \bold{36}.
}

\seealso{
\code{\link{datareplic}}, \code{\link{cirdata}}.
}

\examples{
data(cirgenes)

# These values of kappa estimations were obtained by using the calculations
# explained in Fern?ndez et al (2012).

kappas <- c(2.64773, 3.24742, 2.15936, 4.15314, 4.54357,
               29.07610, 6.51408, 14.19445, 5.66920, 11.12889)

allresults <- list()
resultIsoCIRE <- matrix(ncol = ncol(cirgenes), nrow = nrow(cirgenes))

SCEs <- vector(mode = "numeric", length = nrow(cirgenes))
pvalues <- vector(mode = "numeric", length = nrow(cirgenes))

for (i in 1 : nrow(cirgenes)) {
    k <- kappas[i]
    genes <- as.numeric(cirgenes[i, ])
    allresults[[i]] <- cond.test(genes, groups=c(1:length(genes)), kappa = k)
    resultIsoCIRE[i, !is.na(cirgenes[i, ])] <- round(unlist(allresults[[i]]$CIRE), digits = 3)
    SCEs[i] <- allresults[[i]]$SCE
    pvalues[i] <- allresults[[i]]$pvalue
 }

# resultIsoCIRE is a matrix with the CIRE for each experiment
# pvalues is a vector with the result of the test

}

\keyword{datasets}
\keyword{circular}
\keyword{genes}