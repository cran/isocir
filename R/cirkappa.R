"cirkappa"<-function(data, biascorr=TRUE){
# DATE WRITTEN: 4 Ago 2010          LAST REVISED:  12 Ene 2010
# AUTHOR: Sandra Barragan based on the SAS routines written by Miguel A. Fernandez.
# DESCRIPTION: It calculates the estimation of the kappa parameter of a von Mises.
# data should be a matrix (individuals x replications)
# it needs the package circular since the function A1inv
# REFERENCE: Mardia, K. and Jupp, P. (2000). Directional Statistics.
# SEE ALSO: mrl, CTi, CIREi.
  a<-apply(data,1,mean,na.rm=TRUE)
  data<-data[!is.na(a),]
  b<-apply(data,2,mean,na.rm=TRUE)
  data<-data[,!is.na(b)]
  n<-ncol(data)*nrow(data)
  A<-(1/n)*sum(sin(data),na.rm=TRUE)
  B<-(1/n)*sum(cos(data),na.rm=TRUE)
  R<-sqrt((A^2)+(B^2))
  k<-A1inv(R)
  if(biascorr == TRUE) {
    if(k < 2){
      k<-max(k-(2/(n*k)), 0)
      }else{
        k<-(((n-1)^3)*k)/((n^3)+n)
        }
      }
  return(k)
}