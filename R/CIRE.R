"CIRE"<-function(data, groups=NULL, circular=TRUE){

# DATE WRITTEN: 04 Mar 2010          LAST REVISED:  07 Feb 2013
# AUTHOR: Sandra Barragan based on the SAS routines written by Miguel A. Fernandez.
# DESCRIPTION: Computes the algorithm to obtain the CIRE.
# REFERENCE: Rueda, C., Fernandez, M. A. and Shyamal, P. (2009). 
#             Estimation of parameters subject to order restrictions 
#             on a circle with application to estimation of phase angles of cell-cycle genes. 
#             JASA.
# SEE ALSO: cond.test, isocir.

  
  

##################################################
##  auxiliary functions needed for CIRE:
##################################################

rotation <- function(ordenation, CIRE){
  # rota el orden introducido en el argumento "ordenation" y CIRE
  # hasta que el primer elemento sea el 1
  while(ordenation[1] != 1){
    ordenation <- c(ordenation[2:length(ordenation)], ordenation[1])
    CIRE <- c(CIRE[2:length(CIRE)], CIRE[1])
  }
  return(list(ordenation=ordenation, CIRE=CIRE))
}

 
  cirmean<-function (data) {
     data <- data[!is.na(data)]
     A <- sum(sin(data))/length(data)
    B <- sum(cos(data))/length(data)
if ((A >= 0) && (B > 0)) {
  result <- atan(A/B)
}
if((A>0) && (B==0)){result <- pi/2}
if (B < 0) {
  result <- atan(A/B) + pi
}
if ((A < 0) && (B >= 0)) {
  result <- atan(A/B) + (2 * pi)
}
     return(result)
 }
 


mmeans <- function(data){
   means <- matrix(0, nrow=length(data), ncol=length(data))
   for(l in 1:length(data)){
    for(m in l:length(data)){
     means[l,m]<-cirmean(data[l:m])
    }
   }
   means<-means+t(means)
   diagonal<-diag(diag(means),nrow=nrow(means),ncol=ncol(means))
   means<-means-(diagonal/2)
   return(means)
   # the result is a matrix nxn where n=length(data)
}

### the main function:

funCIRE <- function(data, means, mrl, circular){
  q <- length(data)
  if(circular){
     fis.tilde<-.Call("funCIREnuevo", data, means, mrl, PACKAGE="isocir")
       CIRE<-fis.tilde[1:q]
       SCE<-fis.tilde[q+1]   
   } #fin circular
   
   if(!circular){
     data1 <- data
    
     lA1.PAVA <- list()
     lA2.PAVA <- list()
     lA3.PAVA <- list()
     
     resultlist<-.Call("funCIRE", data, means, mrl, PACKAGE="isocir")
     
     lA1.PAVA<-resultlist[[1]]
     lA2.PAVA<-resultlist[[2]]
     lA3.PAVA<-resultlist[[3]]
     fis.tilde<-matrix(ncol=length(data), nrow=(length(data)+1)*length(data))
     z<-1
   
   fis.tilde[z,]<-lA3.PAVA[[1]]
   z<-z+1
   for(j in 2:length(data)){
     fis.tilde[z,]<-c(lA2.PAVA[[1]][[j-1]],lA3.PAVA[[j]])
     z<-z+1     
   }
   fis.tilde[z,]<-lA2.PAVA[[1]][[length(data)]]
   z<-z+1
   for(i in 2:length(data)){
     for(j in i:length(data)){
       if(j!=length(data)){
         fis.tilde[z,]<-c(lA1.PAVA[[i-1]],lA2.PAVA[[i]][[j]],lA3.PAVA[[j+1]])
         z<-z+1       
       }
       if(j==length(data)){
         fis.tilde[z,]<-c(lA1.PAVA[[i-1]],lA2.PAVA[[i]][[j]])
         z<-z+1         
       }
     }     
   }
   for(i in 1:(length(data)-1)){
     fis.tilde[z,]<-c(lA1.PAVA[[i]],lA3.PAVA[[i+1]])
     z<-z+1
   }
   fis.tilde[z,]<-c(lA1.PAVA[[length(data)]])
   fis.tilde<-fis.tilde[1:z,]
   addition<-apply(fis.tilde,1,sum)
   if(all(!addition<10000)){
     SCE <- 10000
     CIRE <- rep(10000,q)
   }
   if(!all(!addition<10000)){
     fis.tilde<-rbind(fis.tilde[addition<10000,])
     SCE<-apply(fis.tilde,1,sce,data,mrl)
     CIRE<-fis.tilde[which.min(SCE),]
     SCE<-min(SCE)
   }
   rm(fis.tilde,lA1.PAVA,lA2.PAVA,lA3.PAVA)
     
     
   } # find no circular
     
     
return(list(CIRE=CIRE,SCE=SCE))
}# end function funCIRE

####################

##################################################
#####  control of arguments: ####
##################################################

if(missing(data)){stop("Data cannot be misssing")}
if(all(c(!is.matrix(data),!is.vector(data),
         !(is.circular(data)&&is.numeric(data)),!is.data.frame(data)))){
      stop("The argument data must be a matrix or a vector")}
if(length(data) == 1){stop("Data cannot be of length 1")}
#if(!is.matrix(data)){ if(is.vector(data)){data<-cbind(data)} }
if(!is.numeric(data)||!is.matrix(data)) {
  if (is.vector(data)) {
    data <- cbind(as.numeric(data))
  }
  if (is.data.frame(data)) {
    if(nrow(data)>1){data <- as.matrix(data)}
    if(nrow(data)==1){data <- cbind(as.numeric(data))}
  }
}


means <- NULL
for(i in 1:nrow(data)){
  if(all(is.na(data[i,]))){
    means[i] <- NA
  }
  else{
    means[i] <- cirmean(data[i,])
  }  
}
  
#apply(data,1,cirmean)
if(is.null(groups)){groups <- c(1:length(means))}
resultant <- mrl(data)[!is.na(means)]
#groups <- order(order(groups[which(!is.na(means))]))
groups<- groups[!is.na(means)]
means <- means[!is.na(means)]

### ajustamos los datos para que no haya 0s o 2pis exactos

for(i in 1:length(means)){
  if(means[i]==0){means[i]<-1e-50}
  if(means[i]==(2*pi)){means[i]<-6.28318}
}


#comprobamos que al quitar los NAs no se haya quedado sin algun grupo
comprobar<-c(1:nrow(data))
donde<-rep(NA,length(comprobar))
for(i in 1:length(comprobar)){
  if(length(which(groups==comprobar[i]))>0){donde[i] <- i}
}

if(any(is.na(donde))){ # si se ha quedado sin algun grupo:
  falta<-sort(which(is.na(donde)))
  groups[which(groups>falta[1])]<-groups[which(groups>falta[1])]-1  
  falta<-falta-1
  if(length(falta)>1){ 
    for(i in 2:length(falta)){
      groups[which(groups>falta[i])]<-groups[which(groups>falta[i])]-1
      falta<-falta-1
    }
  }
}


n <- as.vector(table(groups))
b<-as.numeric(names(table(groups)))
if(length(b)!=b[length(b)]){stop("There is some missing level")}
k <- length(n)
listdata <- list()
lindex <- list()
length(listdata) <- length(lindex) <- k
for(i in 1:length(groups)){
  if(!is.null(listdata[[groups[i]]])){
    listdata[[groups[i]]] <- c(listdata[[groups[i]]],means[i])
    lindex[[groups[i]]] <- c(lindex[[groups[i]]], i)  
  }
  if(is.null(listdata[[groups[i]]])){
    listdata[[groups[i]]] <- means[i]
    lindex[[groups[i]]] <- i
  }
}
# ordresult <- groups
# names(ordresult) <- resultant
# sortresult <- as.numeric(names(sort(ordresult)))


# all necessary data:
ldata <- listdata
rm(listdata)
#mrl <- sortresult
num <- n
nr <- prod(factorial(num))
nc <- sum(num)
n <- factorial(num)


# create mdata, the matrix with all the possibilities
# dim>1 just in the case of the partial order
indexes<-matrix(ncol=nc,nrow=nr)
mdata <- matrix(ncol=nc,nrow=nr)
aux <- 1
for(p in 1:length(ldata)){
 l <- length(ldata[[p]])
 if(l == 1){
  mdata[, aux] <- rep(ldata[[p]], nr)
  indexes[, aux] <- rep(aux, nr)
  aux <- aux + 1 
  }# end if group of size 1
 if(l > 1){
  iaux <- c(aux:(aux + l - 1))
  permu <- permn(ldata[[p]])
  ipermu <- permn(iaux)
  lp <- length(permu)
  nrep <- nr / lp
  if(p == 1){
   maux <- matrix(nrow = lp, ncol = l)
   miaux <- matrix(nrow = lp, ncol = l)
   for(d in 1:lp){
    maux[d,] <- permu[[d]]
    miaux[d,] <- ipermu[[d]]
    }
   for(c in 1:nrep){
    mdata[(n[1]*(c - 1) + 1):(c*n[1]), 1:l] <- maux
    indexes[(n[1]*(c - 1) + 1):(c*n[1]), 1:l] <- miaux
    }
   aux <- aux+l
   } # close if p=1
  if(p != 1){
   repe <- prod(n[1:(p - 1)])
   maux <- matrix(nrow = repe, ncol = l)
   miaux <- matrix(nrow = repe, ncol = l)
   for(a in 1:(nr / repe)){
    aux1 <- a%%n[p]
    if(aux1 == 0){aux1 <- n[p]}
    repe <- prod(n[1:(p - 1)])
    maux <- matrix(nrow = repe, ncol = l)
    miaux <- matrix(nrow = repe, ncol = l)
    for(d in 1:repe){
     maux[d,] <- permu[[aux1]]
     miaux[d,] <- ipermu[[aux1]]
     }
    if(a > n[p]){aux1 <- aux1 + 1}
    mdata[(((a - 1)*repe) + 1):(a*repe), aux:(aux + l - 1)] <- maux
    indexes[(((a - 1)*repe) + 1):(a*repe), aux:(aux + l - 1)] <- miaux
    } # close for a
   aux <- aux + l
   } # close if p!=1
  } # close if l>1
 } # close p

# coger los mrl en el orden correcto:
mrls<-matrix(ncol = ncol(mdata), nrow = nrow(mdata))
for(i in 1:nrow(mdata)){mrls[i,]<-resultant[indexes[i,]]}


CIRE <- matrix(ncol = ncol(mdata), nrow = nrow(mdata))
SCE <- array()

#### calculate the CIRE ##### 
# for each row of mdata
for(m in 1:nrow(mdata)){
 data <- mdata[m,]
 mrl2 <- mrls[m,]
 ind <- indexes[m,]
 means <- mmeans(data)
 resultfunCIRE <- funCIRE(data, means, mrl2, circular)
 CIRE[m,] <- resultfunCIRE$CIRE
 SCE[m] <- resultfunCIRE$SCE
 rm(resultfunCIRE)
 # In case of circular order:
 if(circular == TRUE){
  aux <- 1
  for (i in 2:length(data)){
   data <- c(data[2:length(data)], data[1])
   mrl2 <- c(mrl2[2:length(mrl2)], mrl2[1])
   ind <- c(ind[2:length(ind)], ind[1])
   means <- mmeans(data)
   resultfunCIRE <- funCIRE(data, means, mrl2, circular)
   # If the SCE is better, we take that new result
   if(resultfunCIRE$SCE < SCE[m]){
    CIRE[m,] <- resultfunCIRE$CIRE
    SCE[m] <- resultfunCIRE$SCE
    rm(resultfunCIRE)
    indexes[m,] <- ind
    aux <- i
    } # close if<SCE
   } # close for(i)
  } # close circular=TRUE
 }# close for(m)

rm(mdata, means)

CIREfinal <- unique(rbind(CIRE[c(1:length(SCE))[SCE == min(SCE)],]),MARGIN=1)
minSCE <- unique(SCE[SCE == min(SCE)])
finalindexes <- rbind(indexes[c(1:length(SCE))[SCE == min(SCE)],])[1,]

aux <- rotation(ordenation=finalindexes, CIRE=CIREfinal)
rotCIRE<-aux$CIRE[aux$ordenation]

iniorder <- unlist(lindex)
#finorder <- iniorder[aux$ordenation] # interesante si orden parcial

# devolvemos el CIRE en el orden en el que se han metido los datos
# es decir que pra la misma posicion el valor del CIRE \tilde{theta}_i 
# es el correspondiente al dato \theta_i bajo el orden introducido
# en groups (posiciones del orden)
# no tiene mucho sentido devolver esto y hacerlo como una lista
# ya que eso seria adecuado cuando hubiera ordenes parciales
# y devolverlo en el orden introducido con sus conjuntos.
lCIRE <- list()
# length(lCIRE) <- k
# # no valido en caso de ordenes parciales, se sale de rango
# for(i in 1:length(posorders2)){
#   if(!is.null(lCIRE[[posorders2[i]]])){lCIRE[[posorders2[i]]] <- c(lCIRE[[posorders2[i]]], rotCIRE[i])}
#   if(is.null(lCIRE[[posorders2[i]]])){lCIRE[[posorders2[i]]] <- rotCIRE[i]}
# }

# Para que valga tanto para total como parcial:
finalCIRE <- rotCIRE[order(iniorder)]
length(lCIRE) <- length(finalCIRE)
for(i in 1:length(finalCIRE)){lCIRE[[i]] <- finalCIRE[i]}

# devuelve 

lexit <- list(cirmeans = ldata, SCE = minSCE, CIRE = lCIRE)

class(lexit) <- "isocir"
return(lexit)

} ###### end function CIRE





