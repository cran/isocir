"mrl"<-function(data){
# DATE WRITTEN: 4 Ago 2010          LAST REVISED:  11 Ene 2010
# AUTHOR: Sandra Barragan based on the SAS routines written by Miguel A. Fernandez.
# DESCRIPTION: data is a matrix of data where:
#                the rows are the items and 
#                the colums are the replications.
# REFERENCE: Mardia, K. and Jupp, P. (2000). Directional Statistics.
# SEE ALSO: cirkappa, CTi, CIREi.

A<-array(); B<-array()
for(i in 1:nrow(data)){
 A[i]<-sum(sin(data[i,!is.na(data[i,])]))
 B[i]<-sum(cos(data[i,!is.na(data[i,])])) 
 }
r_g<-sqrt((A^2)+(B^2))
return(r_g)
}