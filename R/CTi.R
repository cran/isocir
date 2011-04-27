"CTi"<-function(data, levels=c(1:nrow(data)), kappa=NULL){
# DATE WRITTEN: 4 Ago 2010          LAST REVISED:  23 Mar 2011
# AUTHOR: Sandra Barragan based on the SAS routines written by Miguel A. Fernandez.
# DESCRIPTION: 
# data are the unordered estimator
# replic=T means that data have replication => kappa estimated
# kappa is the value of kappa when is known
# REFERENCE: Fernandez, M. A., Rueda, C. and Shyamal, P. (2011). 
#             Isotropic order among core set of orthologs conserved between 
#             budding and fission yeasts. Preprint.
# SEE ALSO: cirkappa, mrl, CIREi.

#if(length(levels)==1&&levels[1]==1){levels<-c(1:nrow(data))}
q<-nrow(data)
n<-as.vector(table(levels))
posi<-prod(factorial(n))
resultCIRE<-CIREi(data, levels, isotropic=TRUE)
SCE<-resultCIRE$SCE
ordermeans<-NULL
for(i in 1:length(resultCIRE$CIRE)){
  ordermeans<-c(ordermeans,as.vector(resultCIRE$CIRE[[i]],mode="numeric"))
  }
data1<-sort(ordermeans)
data2<-c(data1[2:q],data1[1])
difference<-ifelse((data1-data2)==0,0,1)
m<-sum(difference)

if(ncol(data)==1){
  kap<-kappa
  T<-2*kap*SCE
  pvalue<-pchisq(T,q-m,lower.tail=FALSE)*(1-((posi)/(factorial(q-1))))
  } # end kappa known

if(ncol(data)>1){
  kap<-cirkappa(data)
  T<-2*kap*SCE
  pvalue<-pf(T,q-m,q-1,lower.tail=FALSE)*(1-((posi)/(factorial(q-1))))
  } # end kappa unknown and replications

return(pvalue)
} # end CTi
