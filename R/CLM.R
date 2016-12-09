"CLM" <- function(data, order0, ws=NULL, control.method=c("msce", "cirktau")){
  
  # DATE WRITTEN: 2 May 2012          LAST REVISED:  7 May 2013
  # AUTHOR: Sandra Barragan.
  # DESCRIPTION: Computes the Circular Local MINIMIZATION (CLM).
  
  cirmean<-function (data) {
    data <- data[!is.na(data)]
    A <- sum(sin(data))
    B <- sum(cos(data))
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
  
  
  cirPAVA <-function(data){
    # DATE WRITTEN: 4 Mar 2010          REWRITTEN:  25 Ene 2013 (totalmente modificada)
    # AUTHOR: Sandra Barragan 
    # DESCRIPTION: It computes PAVA with the circular mean.
    data1<-as.numeric(data)
    problems <- ifelse(data[1:length(data)-1]-data[2:length(data)]>0,1,0)
    while(sum(problems)>0){
      
    for(i in 1:(length(data)-1)){     
      if(problems[i]>0){
        a<- i
        a<- min(c(1:length(data1))[data1==data1[a]])
        b<- i+1
        b<-max(c(1:length(data1))[data1==data1[b]])
        data1[a:b]<-cirmean(data[a:b])        
      }
    }
    problems <- ifelse(data1[1:length(data1)-1]-data1[2:length(data1)]>0,1,0)
    }
    return(data1)
  }
  
  comparison<-function(ordenation, i, j, k, cirorders, data, ws, lmsces, msce){
    
    # hace la comparacion del triplete en las posiciones (i,j,k)
    # del orden ordenation segun los ordenes dados en cirorders
    
    # ordenation is a VECTOR
    exit <- FALSE
    a <- i
    b <- j
    c <- k
    auxorder <- ordenation # the order to be compared 
    # auxorder is a VECTOR y tiene la mejor ordenacion en cada momento 
    best <- msce # esto es para guardar el msce de la ordenacion mejor en cada momento
    bestlist <- lmsces # esto es para guardar la lista de CIREs de la ordenacion mejor en cada momento
    while (exit == FALSE){
      auxorderb <- auxorder #auxorderb sera la ordenacion permutada a comprobar
      auxorderb[c(a,b,c)]<-auxorder[c(a,c,b)] # realizamos la permutacion
      auxlmscesb <-  bestlist # esto es solo porque nos interesa la estructura
      pbest <- 0  # para guardar el msce global de auxorderb (sera un valor)
      pmsces <- NULL # para guardar los msces en cada experimento (sera un vector)
      for(i in 1:nrow(data)){
        CIREi <- as.vector(bestlist[[i]]$CIRE)
        if(any(is.na(c(CIREi[b],CIREi[c])))){
          CIREi[c(b,c)] <- CIREi[c(c,b)]
          auxlmscesb[[i]]$msce <- bestlist[[i]]$msce
          auxlmscesb[[i]]$CIRE <- CIREi
        }
        if(!any(is.na(c(CIREi[b],CIREi[c])))){
          # si no es NA ninguna de la posiciones a intercambiar
          extraction <- sort(unique(c(which(CIREi==CIREi[b]),which(CIREi==CIREi[c]))))
          corte <- FALSE
          q<-1
          while(!corte&&q<length(extraction)){
            r<-extraction[q]+1
            s<-extraction[q+1]-1
            if(r<=s){
              if(!all(is.na(CIREi[r:s]))){
                corte <- TRUE
                extraction <- extraction[c((q+1):length(extraction),1:q)] 
              }
            }
            q <- q+1
          }
        r<-1
          index <- ifelse(((extraction[1]-r)%%length(CIREi))==0,length(CIREi),(extraction[1]-r)%%length(CIREi))
        while(is.na(data[i,auxorder][index])){
          r <- r+1
          index <- ifelse(((extraction[1]-r)%%length(CIREi))==0,length(CIREi),(extraction[1]-r)%%length(CIREi))
        }
        extractionC<-c(index,extraction)
        s <- 1
        index <- ifelse(((extraction[length(extraction)]+s)%%length(CIREi))==0,length(CIREi),(extraction[length(extraction)]+s)%%length(CIREi))
        while(is.na(data[i,auxorder][index])){
          s<-s+1
          index <- ifelse(((extraction[length(extraction)]+s)%%length(CIREi))==0,length(CIREi),(extraction[length(extraction)]+s)%%length(CIREi))
        }  
        extractionC<-c(extractionC,index)
        extractionC<-unique(extractionC)
          
          extraction2 <- extraction
          extraction2[c(which(extraction==b),which(extraction==c))]<-c(c,b)
          paraPAVA <- data[i,auxorder][extraction2]
          
                    
          nozero <- all(order(CIREi[extractionC])==c(1:length(CIREi[extractionC])))
          # (max(paraPAVA)-min(paraPAVA)<pi)&&
          condition <- (max(CIREi[extractionC])-min(CIREi[extractionC])<pi)&&(nozero)
          
          if(condition){         
           CIREi[extraction] <- cirPAVA(paraPAVA)
           if(!is.unsorted(CIREi[c(which.min(CIREi):length(CIREi),1:(which.min(CIREi)-1))], na.rm =TRUE)){
           auxlmscesb[[i]]$msce <- sce(data[i,auxorderb],CIREi)/length(which(!is.na(CIREi)))
           auxlmscesb[[i]]$CIRE <- CIREi
           }
           if(is.unsorted(CIREi[c(which.min(CIREi):length(CIREi),1:(which.min(CIREi)-1))], na.rm =TRUE)){
             condition <- FALSE
           }
         }
        if(!condition){
          aux <- CIRE(data[i,], order(auxorderb))
          MSCE<-(aux$SCE)/length(aux$CIRE)
          auxlmscesb[[i]]$msce <- MSCE
          if(any(is.na(data[i,]))){
            posi <- data[i,auxorderb]
            auxlmscesb[[i]]$CIRE <- rep(NA,ncol(data))
            auxlmscesb[[i]]$CIRE[!is.na(posi)] <- unlist(aux$CIRE)[order(order(auxorderb)[-which(is.na(data[i,]))])]}
          if(!any(is.na(data[i,]))){auxlmscesb[[i]]$CIRE <- unlist(aux$CIRE)[auxorderb]}        
        }
        } # si no es NA ni b ni c         
        
        pmsces <- c(pmsces, auxlmscesb[[i]]$msce)
        pbest <- pbest + ws[i]*auxlmscesb[[i]]$msce
      } # fin de i   

    # Para el triplete dado (abc) hacemos los cambios si necesarios
      
      if(best <= pbest){ # no es necesario cambiar, 
        exit <- TRUE # salimos del bucle y de la funcion
      }     
      if(pbest < best){   # se cambia
        auxorder <- auxorderb # nuestra mejor ordenation ahora es con (acb)
        best <- pbest # este es el mejor msce global
        bestlist <- auxlmscesb # actualizamos la lista de CIREs
        # comprobamos tripletes anteriores si no estamos en el principio
        if(a < (length(ordenation)-1)){ 
          a <- a-1
          b <- b-1
          c <- c-1
        }
        if(a==(length(ordenation)-1)){
          a <- a-1
          b <- b-1
          c <- length(ordenation)
        }
        if(a == length(ordenation)){
          a <- a-1
          b <- length(ordenation)
          c <- 1
        }
        if(a <= 0){ # se ha llegado al principio, finalizamos
          exit <- TRUE
        }
      }      
    } # end while exit FALSE
    return(list(ordenation = auxorder, bestsce = best, lmsces = bestlist))
  } # end comparison SIN EMPATES, atras adelante y + eficiente
  
  taucomparison<-function(ordenation, i, j, k, data, ws, tau){
    
    # hace la comparacion del triplete en las posiciones (i,j,k)
    
    # ordenation is a VECTOR
    exit <- FALSE
    a <- i
    b <- j
    c <- k
    auxorder <- ordenation # the order to be compared #auxorder is a VECTOR
    best <- tau
    while (exit == FALSE){
      auxorderb <- auxorder
      auxorderb[c(a,b,c)]<-auxorder[c(a,c,b)]
      pbest <- mcirktau(data, order(auxorderb), ws)$mtau
      # Para el triplete dado hacemos los cambios si necesarios
      
      if(best >= pbest){  
        # no es necesario cambiar, finalizamos.
        exit <- TRUE
      }     
      if(pbest > best){
        # se cambia
        auxorder <- auxorderb
        best <- pbest
        
        if(a < (length(ordenation)-1)){
          # si vamos hacia atras
          a <- a-1
          b <- b-1
          c <- c-1
        }
        if(a==(length(ordenation)-1)){
          # si vamos hacia atras
          a <- a-1
          b <- b-1
          c <- length(ordenation)
        }
        if(a == length(ordenation)){
          # si vamos hacia atras
          a <- a-1
          b <- length(ordenation)
          c <- 1
        }
        if(a <= 0){ # se ha llegado al principio, finalizamos
          exit <- TRUE
        }
      }
      
      
    } # end while exit FALSE
    return(list(ordenation = auxorder, bestau = best))
  } # end taucomparison
  
  
  
  ##################### MAIN CODE ####
  
  #### control of arguments:
  
  if(!is.matrix(data)&&!is.data.frame(data)){stop("data must be a matrix")}
  if(is.data.frame(data)){data <- as.matrix(data)}
  if(!is.matrix(order0) && !is.vector(order0)){stop("order0 must be a vector or a matrix")}
  if(is.matrix(order0)){order0 <- as.vector(order0)}
  if(length(order0)!=ncol(data)){stop("Error: the dimensions do not match")}
  if(is.null(ws)){ws<-(rep(1,nrow(data)))/nrow(data)}
  if(sum(ws)!=1){ws<-ws/sum(ws)}
  control.method <- match.arg(control.method, c("msce", "cirktau"))
  
  if(control.method=="msce"){
  cirorders <- t(apply(data, 1, order)) # positions in the orders
  posorders <- t(apply(cirorders, 1, order))
  #msce0 <- msce(data, order(order0), ws)$msce
  lmsces <- list(list(msce=0, CIRE=array(NA,dim=ncol(data))))
  length(lmsces) <- nrow(data)
  for(i in 1:(nrow(data)-1)){lmsces[[i+1]] <- list(msce=0, CIRE=array(NA,dim=ncol(data)))}
  msce0 <- 0
  msces2 <- NULL
  for(i in 1:nrow(data)){
    aux <- CIRE(data[i,], order(order0))
    MSCE<-(aux$SCE)/length(aux$CIRE)
    lmsces[[i]]$msce <- MSCE
    if(any(is.na(data[i,]))){
      posi <- data[i,][order0]
      lmsces[[i]]$CIRE[!is.na(posi)] <- unlist(aux$CIRE)[order(order(order0)[-which(is.na(data[i,]))])]
    }
    if(!any(is.na(data[i,]))){lmsces[[i]]$CIRE <- unlist(aux$CIRE)[order0]}
    msces2 <- c(msces2, MSCE)
    msce0 <- msce0 + ws[i]*MSCE
  } 
  
  ne <- ncol(cirorders) # number of elements
  nr <- nrow(cirorders) # number of orders
#   msce1 <- msce(data, order(order0[length(order0):1]), ws)$msce
#   if(msce1 < msce0){
#     msce0 <- msce1
#     order0 <- order0[length(order0):1]
#   }
  ordenation <- order0
  msce <- msce0
  for (l in 1:ne){ # for each element
    i <- l 
    if (i != ne){
      j <- i+1
      k <- i+2}
    if (i == (ne-1)){
      j <- ne
      k <- 1}
    if (i == ne){
      j <- 1
      k <- 2}
    if(l < ne){
      aux <- comparison(ordenation, i, j, k, cirorders, data, ws, lmsces, msce)
      ordenation <- aux$ordenation
      lmsces <- aux$lmsces
      msce <- aux$bestsce
    }
    if(l == ne){
      aux <- comparison(ordenation, i, j, k, cirorders, data, ws, lmsces, msce)
      ordenation <- aux$ordenation
      lmsces <- aux$lmsces
      bestsce <- aux$bestsce
    }
  } # end for l
  
  # devuelve el orden circular fusion de todos los ordenes circulares dados (cirorders)
  return(list(order0 = order0, msce0 = msce0, final_order = ordenation, bestsce = bestsce))
}
##### tau ####
if(control.method=="cirktau"){
  tau0 <- mcirktau(data, order(order0), ws)$mtau
  
  ne <- ncol(data) # number of elements
  nr <- nrow(data) # number of orders
  ordenation <- order0
  tau <- tau0
  for (l in 1:ne){ # for each element
    i <- l 
    if (i != ne){
      j <- i+1
      k <- i+2}
    if (i == (ne-1)){
      j <- ne
      k <- 1}
    if (i == ne){
      j <- 1
      k <- 2}
    if(l < ne){
      aux <- taucomparison(ordenation, i, j, k, data, ws, tau)
      ordenation <- aux$ordenation
      tau <- aux$bestau
    }
    if(l == ne){
      aux <- taucomparison(ordenation, i, j, k, data, ws, tau)
      ordenation <- aux$ordenation
      bestau <- aux$bestau
    }
  } # end for l
  return(list(order0 = order0, tau0 = tau0, final_order = ordenation, bestau = bestau))
  
} #fin tau
  
} # end CLM: Circular Local Minimization
