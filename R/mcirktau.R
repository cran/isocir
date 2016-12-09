"mcirktau" <- function(data, posorder, ws=NULL){
  
  # DATE WRITTEN: 22 Ene 2013          LAST REVISED:  22 Ene 2013
  # AUTHOR: Sandra Barragan.
  # DESCRIPTION: Computes the mean Circular Kendall's tau between
  # a circular order and a data set.
  # REFERENCE:  Fisher (1993) 
  
  ntaus <- apply(data, 1, cirKendall, posorder)
  if(is.null(ws)){ws <- rep(1, nrow(data))}
  return(list(mtau=weighted.mean(ntaus, w=ws, na.rm=TRUE), ntaus))
}