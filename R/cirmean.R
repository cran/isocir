"cirmean"<-function(data){
# DATE WRITTEN: 4 Mar 2010          LAST REVISED:  11 Ene 2010
# AUTHOR: Sandra Barragan based on the SAS routines written by Miguel A. Fernandez.
# DESCRIPTION: Computes the circular mean.
# REFERENCE: Mardia, K. and Jupp, P. (2000). Directional Statistics.
# SEE ALSO: CIREi, cirPAVA, cirSCE.
 
 data<-data[!is.na(data)]
 A<-sum(sin(data))
 B<-sum(cos(data)) 
 if((A>=0)&&(B>0)){result<-atan(A/B)}
 if(B<0){result<-atan(A/B)+pi}
 if((A<0)&&(B>0)){result<-atan(A/B)+(2*pi)}
 return(result)
 } # close cirmean