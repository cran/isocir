"mrl"<-function(data){
  # DATE WRITTEN: 4 Ago 2010          LAST REVISED:  22 Ene 2013
  # AUTHOR: Sandra Barragan.
  # DESCRIPTION: data is a matrix of data where:
  #                the rows are the items and 
  #                the colums are the replications.
  # REFERENCE: Mardia, K. and Jupp, P. (2000). Directional Statistics.
  # SEE ALSO: cirkappa, cond.test, CIRE.
  
  # control of arguments:
  if(!is.matrix(data)&&!is.data.frame(data)){
    if(is.vector(data)){data<-cbind(data)}
    else{stop("data must be a matrix or a vector")}
  }
  
  ####
  A<-array(); B<-array()
  for(i in 1:nrow(data)){
    A[i]<-sum(sin(data[i,]),na.rm=TRUE)
    B[i]<-sum(cos(data[i,]),na.rm=TRUE) 
  }
  r_g<-sqrt((A^2)+(B^2))/ncol(data)
  return(r_g)
  
# by using the function rho.circular from circular package:
# if(all(c(!is.matrix(data),!is.vector(data),!is.data.frame(data)))){stop("data must be a matrix or a vector")}
# if(!is.circular(data)){data <- suppressWarnings(as.circular(data))}
# 
# if(is.matrix(data)){ r <- suppressWarnings(apply(data,1,rho.circular))}
# if(is.null(dim(data))){r <- suppressWarnings(rho.circular(data))}
# return(r)
  
}