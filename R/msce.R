

"msce" <- function(data, posorder, ws=NULL, ...){
  # DATE WRITTEN: 22 Ene 2013          LAST REVISED:  22 Ene 2013
  # AUTHOR: Sandra Barragan.
  # DESCRIPTION: Computes the MSCE between and order and a set of data.
  # REFERENCE: 
  
  msces <- NULL
  for(j in 1:nrow(data)){
    msces[j] <- (CIRE(data[j,], groups=posorder, ...)$SCE)/length(data[j,!is.na(data[j,])])
  }
  if(is.null(ws)){ws <- (rep(1, nrow(data)))/nrow(data)}
  return(list(msce=weighted.mean(msces, w=ws),msces))  
}