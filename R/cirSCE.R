"cirSCE"<-function(arg1,arg2,mrl=1){
# DATE WRITTEN: 4 Mar 2010          LAST REVISED:  12 Abr 2011
# AUTHOR: Sandra Barragan based on the SAS routines written by Miguel A. Fernandez.
# DESCRIPTION: It computes the Sum of Circular Errors of two vectors
# REFERENCE: Mardia, K. and Jupp, P. (2000). Directional Statistics.
# SEE ALSO: CIREi, cirPAVA, cirmean.
if(is.vector(arg2)){
 proveNA1<-rep(0,length(arg1))
 proveNA1[is.na(arg1)]<-1
 if(sum(proveNA1)>=1){
  arg1<-arg1[proveNA1==0]
  arg2<-arg2[proveNA1==0]
  }
 proveNA2<-rep(0,length(arg2))
 proveNA2[is.na(arg2)]<-1 
 if(sum(proveNA2)>=1){
  arg2<-arg2[proveNA2==0]
  arg1<-arg1[proveNA2==0]
  } 
 point1<-arg1
 point2<-arg2
 if(length(mrl)>1){mrls<-mrl}
 if(length(mrl)==1){mrls<-rep(1,length(point1))}
 } # end if arg2 vector
if(is.matrix(arg2)){
 point1<-arg1
 point2<-apply(arg2,1,cirmean) 
 if(length(mrl)>1){mrls<-mrl}
 if(length(mrl)==1){mrls<-mrl(arg2)}
 } # end if arg2 matrix
aux<-mrl*(1-cos(point1-point2))
SCE<-sum(aux)
return(SCE)
}