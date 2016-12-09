"cond.test"<-function(data, groups=NULL, kappa = NULL, biasCorrect = TRUE){
# DATE WRITTEN: 4 Ago 2010          LAST REVISED:  22 Ene 2013
# AUTHOR: Sandra Barragan based on the SAS routines written by Miguel A. Fernandez.
# DESCRIPTION: 
# data are the unordered estimator
# replic=T means that data have replication => kappa estimated
# kappa is the value of kappa when is known
# REFERENCE: Fernandez, M. A., Rueda, C. and Shyamal, P. (2012). 
#              Identification of a core set of signature cell cycle genes whose relative order of time to peak expression is conserved across
#              species, \emph{Nucl. Acids Res.} \bold{40}, n7: pp 2823--2832. doi:10.1093/nar/gkr1077.
# SEE ALSO: mrl, CIRE.

if(is.vector(data)){data <- cbind(data)}
if(ncol(data) == 1 && missing(kappa)){stop("If data do not have replications, kappa must be known")}
if(is.null(groups)){stop("The order in the null hypothesis must be defined in the argument groups")}

n <- as.vector(table(groups[!apply(data,1,is.na)]))
posi <- prod(factorial(n))
resultCIRE <- CIRE(data, groups, circular=TRUE)
SCE <- resultCIRE$SCE
ordermeans <- unlist(resultCIRE$CIRE)
data1 <- sort(ordermeans)
data2 <- c(data1[2:length(data1)],data1[1])
difference <- ifelse((data1 - data2) == 0, 0, 1)
m <- sum(difference)
q <- length(ordermeans)

if(!is.null(kappa)& SCE!=0){
  kap <- kappa
  ST <- 2*kap*SCE
  pvalue <- pchisq(ST, q-m, lower.tail=FALSE)*(1 - ((posi)/(factorial(q - 1))))
  } # end kappa known

if(is.null(kappa)){
  Rs<-apply(data,1,mrl)
  kappa<-A1inv(mean(Rs))
  kap <- kappa
  if(biasCorrect == TRUE) {
    N<-ncol(data)*nrow(data)
    if(kap < 2){
      kap<-max(kap-(2/(N*kap)), 0)
    }else{
      kap<-(((N-1)^3)*kap)/((N^3)+N)
    }
  }
if(SCE!=0){
  ST <- 2*kap*SCE
  pvalue <- pf(ST, q-m, q-1, lower.tail=FALSE)*(1 - ((posi)/(factorial(q - 1))))
}
} # end kappa unknown and replications
if(SCE==0){
  pvalue<-1
  #warning("SCE=0")
}
lexit <- resultCIRE
class(lexit)<-"isocir"
lexit$pvalue <- pvalue
lexit$kappa <- kappa
attr(lexit,"estkappa")<-ifelse(is.null(kappa),"Kappa has been estimated","Kappa was known")

return(lexit)
} # end cond.test
