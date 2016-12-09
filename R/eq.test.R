"eq.test"<-function(data, popu, ws=NULL, method=NULL, control.method=NULL, output=NULL, coef=1, N=500){
  
  # DATE WRITTEN: May 2013          LAST REVISED:  22 Mar 2013
  # AUTHOR: Sandra Barragan.
  # DESCRIPTION: Computes the test of equality of circular orders.
  
  #   data: matrix with the experiments in the row and the elements in the columns.
  #   popu: vector of same length as rows in data, has the number of the population each experiment belongs.
  #   
  
  
  ##################################################
  ##  auxiliary functions needed for eq.test:
  ##################################################
  
  
  rotation <- function(ordenation){
    # rota el orden introducido en el argumento "ordenation" 
    # hasta que el primer elemento sea el 1
    while(ordenation[1] != 1){
      ordenation <- c(ordenation[2:length(ordenation)], ordenation[1])
    }
    return(ordenation)
  }
  
  cirmean<-function (data, ws=NULL) {
    # data is a vector
    if(is.null(ws)){ws<-rep(1,length(data))}
    if(any(is.na(data))){
      ws<-ws[complete.cases(data)]
      ws <- ws/sum(ws)}
    data <- data[complete.cases(data)]
    A <- sum(ws*sin(data))
    B <- sum(ws*cos(data))
    if ((A >= 0) && (B >= 0)) {
      result <- atan(A/B)
    }
    if (B < 0) {
      result <- atan(A/B) + pi
    }
    if ((A < 0) && (B > 0)) {
      result <- atan(A/B) + (2 * pi)
    }
    return(result)
  }
  
  
  ##################################################
  ##  control of arguments ####
  ##################################################
  
  ##data
  if(!is.matrix(data)&&!is.data.frame(data)){stop("Data must be a matrix or data.frame")}
  
  ## ws
  if(is.null(ws)){ws <- rep(1,nrow(data))}
  if(length(ws)!=nrow(data)){stop("Do not match the experiment weights")}
  ## popu
  if(length(popu)!=nrow(data)){stop("Do not match the population info")}
  
  #coef
  if(is.null(coef)){coef<-1}
  
  #method
  if(is.null(method)){method<-"TSP"}
  
  #control.method
  if(is.null(control.method)){control.method<-"alpha3"}
  
  ##################################################
  ##  main code:
  ##################################################
  
  wsg<-ws/sum(ws)
  odata <- NULL
  ows <- NULL
  S <- length(unique(popu))
  
  datos <- list()
  weight <- list()
  auxweight <- NULL
  MSCE <- rep(NA,S)
  probabilities <- NULL
  
  global1 <- ACO(data, method, control.method, wsg, coef)
  j<-0
  while(j<5){
    global1b <- ACO(data, method, control.method, wsg, coef)
    if(global1b$msce<global1$msce){global1 <- global1b}
    j<-j+1
  }
  global2 <- CLM(data, order0=global1$aggre_order, wsg)
  MSCEG <- global2$bestsce
  
  for(i in 1:S){
    odata <- rbind(odata, as.matrix(data[which(popu==i),]))
    datos[[i]] <- rbind(data[which(popu==i),])
    wsi <- ws[which(popu==i)]
    weight[[i]] <- wsi/sum(wsi)
    auxweight[i] <- sum(wsi)
    msceglobal <- msce(datos[[i]], order(global2$final_order), weight[[i]])$msce
    part1 <- ACO(datos[[i]], method, control.method, ws=weight[[i]], coef)
    j<-0
    while(j<5){
      part1b <- ACO(datos[[i]], method, control.method, ws=weight[[i]], coef)
      if(part1b$msce<part1$msce){part1 <- part1b}
      j <- j+1
    }
    if(part1$msce<=msceglobal){order0 <- part1$aggre_order} # Option 2
    if(part1$msce>msceglobal){order0 <- global2$final_order}
    part2 <- CLM(datos[[i]], order0=order0, ws=weight[[i]])
    MSCE[i] <- part2$bestsce
    #if(auxmsce<MSCE[i]){MSCE[i]<-auxmsce}
    probabilities <- c(probabilities, rep((1/(S*length(which(popu==i)))),times=length(which(popu==i))))
  }
  
  
  Tob <- (MSCEG-weighted.mean(MSCE,auxweight))/MSCEG
  globals <- rbind(rotation(global2$final_order))
  opopu <- sort(popu)
  ST <- rep(NA,N)
  
  for(i in 1:N){
    cat(i)
    cat(",")
    MSCEi <- NULL
    permutation <- sample(nrow(odata), size=nrow(odata), replace=TRUE, prob=probabilities)
    datai <- odata[permutation,]
    wsi <- ws[permutation]
    aux <- CLM(datai, order0=ACO(datai, method, control.method, ws=wsi/sum(wsi), coef)$aggre_order, ws=wsi/sum(wsi))
    pesosi <- NULL
    for(j in 1:S){
      auxws <- wsi[which(opopu==j)]
      pesosi[j] <- sum(auxws)
      wsij <- auxws/sum(auxws)
      msceglobal <- msce(datai[which(opopu==j),], order(aux$final_order),wsij)$msce
      orderest1 <- ACO(datai[which(opopu==j),], method, control.method, ws=wsij, coef)
      if(orderest1$msce<=msceglobal){order0 <- orderest1$aggre_order} # Option 2
      if(orderest1$msce>msceglobal){order0 <- aux$final_order}
      MSCEi[j] <- CLM(datai[which(opopu==j),], order0=order0, ws=wsij)$bestsce
      #if(auxmsce<MSCEi[j]){MSCEi[j] <- auxmsce} # Option 1
    }
    
    ST[i] <- ifelse(aux$bestsce==0,100,(aux$bestsce-weighted.mean(MSCEi,pesosi))/aux$bestsce)
    globals <- rbind(globals, rotation(aux$final_order))
  }
  globalsjoint <- cbind(paste(globals[,1],globals[,2],globals[,3],sep=""))
  for(i in 4:ncol(data)){
    globalsjoint <- cbind(paste(globalsjoint[,1],globals[,i],sep=""))
  }
  MFO <- names(which.max(table(globalsjoint)))
  CC <- CCMFO <- (max(table(globalsjoint))/(N+1))*100
  neworder <- globals[which(globalsjoint==MFO)[1],]
  if(MFO!=globalsjoint[1,]){    
    mscenew <- msce(data, order(neworder), ws)$msce
    if(mscenew>=MSCEG){
      CC <- ((table(globalsjoint)[which(names(table(globalsjoint))==globalsjoint[1,])])/(N+1))*100
    }
    if(mscenew<MSCEG){
      globals[1,] <- neworder
      MSCEG <- msce(data, order(neworder), ws)$msce
      Tob <- (MSCEG-weighted.mean(MSCE,auxweight))/MSCEG
    }    
  }
  
  
  cat("\n pvalue=")
  pvalue<-length(which(ST>Tob))/length(which(!is.na(ST)))
  cat(pvalue)
  cat("\n")
  global_order <- globals[1,]
  globals <- cbind(c(Tob,ST),globals)
  
  if(!is.null(output)){
    write.csv2(globals,paste(output,"globalorders.csv",sep=""))
    write.csv2(table(globalsjoint),paste(output,"frequencydist.csv",sep=""))
  }
  
  return(list(allorders=globals, global_order=global_order, pvalue=pvalue, CC=CC, MFO=neworder, CCMFO=CCMFO))
  
} # end eq.test
