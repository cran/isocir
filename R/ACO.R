
"ACO" <- function(data, method=c("Naive", "CB", "CMC", "TSP", "CH"), control.method=c("tau", "MSCE", "pos", "cirmean", "cirmed", "1", "2", "3", "4m", "4c", "bin", "pos", "alpha1", "alpha2", "alpha3", "alpha4", "alphainf", "time", "arc", "chord", "bin", "pos", "cos", "cmean", "mrl", "e3", "ave", "qua", "nat", "natp", "natb"), ws=NULL, coef=1){
  
  # DATE WRITTEN: 20 Sep 2013          LAST REVISED:  20 Nov 2014
  # AUTHOR: Sandra Barragan.
  # DESCRIPTION: Computes the Aggregation Circular Order.
  # REFERENCE: 
 
  
  ###################### auxiliary functions ###################### 
 
  circularOrderFussion <- function(data, type = c(0,1,2,3,4), control.checking = c("count", "majority", "MSCE", "tau"), ws=NULL){
   
    if(type!=4 && type!=3 && type!=2 && type!=1 && type!=0){stop("type must be 0, 1, 2, 3 or 4.")}
    
    cirorders <- t(apply(data,1,order))   
    nr <- nrow(cirorders)
    ne <- ncol(cirorders)
    
    if(type==1 || type==2 || type==3 ||type==4){ # en caso de haber escogido uno de los enfoques con cadenas de Markov
      M <- matrix(ncol=ne,nrow=ne) # transition matrix
      
      ## MC 1 (adaptado al circulo)
      # existe probabilidad positiva e igual en todos los casos
      # (incluso en el de permanecer) de pasar de i a j para todos 
      # aquellos elementos j que en ALGUN orden sea mas corto 
      # el camino de i a j que de j a i manteniendo la direccion del c?rculo.
      if(type == 1){ 
        for (i in 1:ne){
          pre <- NULL
          for (k in 1:nr){
            auxorder <- cirorders[k,]
            pre <- c(pre, auxorder[1:which(auxorder==i)])
          }
          pre <- unique(pre)
          prob <- 1/length(pre)
          M[i,pre] <- prob
          M[i,which(is.na(M[i,]))] <- 0
        } # end for
      } # end type 1
      
      ## MC 2
      if(type == 2){
        div <- floor((ne-1)/2)*nr
        for (i in 1:(ne-1)){
          for (j in (i+1):ne){
            dif <- NULL
            for (k in 1:nr){
              auxorder <- cirorders[k,]
              pos1 <- which(auxorder == i)
              pos2 <- which(auxorder == j)
              if (pos2 > pos1){dif[k] <- (pos2-pos1)-1}
              if (pos2 < pos1){dif[k] <- ((ne-pos1)+pos2)-1}
            }
            dif[which(dif <= (ne-2)/2)] <- 0
            dif[which(dif > (ne-2)/2)] <- 1
            M[i,j] <- sum(dif)/div
            M[j,i] <- length(which(dif==0))/div
          } #end for j
          M[i,i] <- 1 - sum(M[i,], na.rm = TRUE)
        } # end for i
        M[ne,ne] <- 1 - sum(M[ne,], na.rm = TRUE)
      } # end type 2
      
      ## MC 3
      if(type == 3){
        div <- ne*nr
        for (i in 1:(ne-1)){
          for (j in (i+1):ne){
            dif <- NULL
            for (k in 1:nr){
              auxorder <- cirorders[k,]
              pos1 <- which(auxorder == i)
              pos2 <- which(auxorder == j)
              if (pos2 > pos1){dif[k] <- (pos2-pos1)-1}
              if (pos2 < pos1){dif[k] <- ((ne-pos1)+pos2)-1}
            }
            dif[which(dif <= (ne-2)/2)] <- 0
            dif[which(dif > (ne-2)/2)] <- 1
            M[i,j] <- sum(dif)/div
            M[j,i] <- length(which(dif==0))/div
          } #end for j
          M[i,i] <- 1 - sum(M[i,], na.rm = TRUE)
        } # end for i
        M[ne,ne] <- 1 - sum(M[ne,], na.rm = TRUE)
      } # end type 3
      
      ## MC 4 (adaptado al circulo)
      # Existe probabilidad positiva de pasar de i a j 
      # para todos los elementos j que en LA MAYORIA 
      # de los ordenes sea mas corto el camino de i a j que de j a i
      # manteniendo la direccion del c?rculo.
      if(type == 4){
        control.checking <- match.arg(control.checking)
        for (i in 1:ne){
          if(i!=ne){
            for (j in (i+1):ne){
              dif <- NULL
              for (k in 1:nr){
                auxorder <- cirorders[k,]
                pos1 <- which(auxorder == i)
                pos2 <- which(auxorder == j)
                if (pos2 > pos1){dif[k] <- (pos2-pos1)-1}
                if (pos2 < pos1){dif[k] <- ((ne-pos1)+pos2)-1}
              }
              if(control.checking == "majority"){
                dif[which(dif <= (ne-2)/2)] <- 0
                dif[which(dif > (ne-2)/2)] <- 1
                eq <- nr/2
              }
              if(control.checking == "count"){eq <- ((ne-2)*nr)/2}
              if (sum(dif) <= eq){
                M[i,j] <- 0
                M[j,i] <- 1/ne
              }
              if (sum(dif) > eq){
                M[i,j] <- 1/ne
                M[j,i] <- 0
              }    			
            } # end for j
          } # end if i!=ne
          M[i,i] <- 1 - sum(M[i,],na.rm = TRUE)		
        } # end for i 
      } # end if type 4
      
      descomp <- eigen(t(M), symmetric = FALSE)
      eigvect <- Re(descomp$vectors[,1])
      if(all(eigvect<0)){eigvect <- abs(eigvect)}
      aggre_order <- order(eigvect, decreasing=TRUE)
      tau <- mcirktau(data, order(aggre_order), ws)$mtau 
      sce <- msce(data,  order(aggre_order), ws)$msce 
    } # end MC approaches
    
    ## Metodo Naive
    # Se escoge aquel orden de entre los dados cuya/o
    # Tau Circular de Kendall media/ MSCE sea mayor respecto 
    # al resto de los ordenes circulares dados
    if(type==0){
      if(control.checking=="MSCE"){
      scenaive <- rep(NA, nrow(data))
      for(i in 1:nrow(data)){
        scenaive[i]<-msce(data, order(order(data[i,])), ws)$msce        
      }
      sce<-min(scenaive, na.rm=TRUE)
      aggre_order<-order(data[which.min(scenaive),])
      tau <- mcirktau(data, order(aggre_order), ws)$mtau
      }
      
      if(control.checking=="tau"){
        taus <- rep(NA, nrow(data))
      for (i in 1:nrow(data)){
       taus[i] <- mcirktau(data, order(order(data[i,])), ws)$mtau
      }
      tau <- max(taus)
      aggre_order<- order(data[which.max(taus),])
      sce <- msce(data, order(aggre_order), ws)$msce
      }
    }
        
    return(list(aggre_order = aggre_order, msce=sce, mtau = tau))
  } # end circularOrderFussion
  
  
  CORAM <- function(data, ws, option = c("alpha3","bin", "pos", "alpha1", "alpha2",  "alpha4", "alphainf", "time", "arc", "chord", "dR-dL")){
    
    option <- match.arg(option)  
    cirorders <- matrix(ncol=ncol(data),nrow=nrow(data))
    for(i in 1:nrow(data)){
      cirorders[i,!is.na(data[i,])] <- order(data[i,!is.na(data[i,])])
    }
#    cirorders <- t(apply(data, 1, order))
    
    if(option == "arc"){
      
      ####  Con el arco unidireccionalmente
      # -------------------------
      
      mat <- list()
      n <- ncol(data)
      for (i in 1:nrow(data)){
        transition <- array(NA, dim=c(n,n))
        
        for(j in 1:n){
          if(j != n){
            for(l in (j+1):n){
              if(!is.na(data[i,j])&&!is.na(data[i,l])){
                aux1 <- (data[i,l] - data[i,j])%%(2*pi)
                aux2 <- (data[i,j] - data[i,l])%%(2*pi)
                transition[j,l] <- aux1
                transition[l,j] <- aux2 
              }
            } # recorrer cada fila
          } # if j!=n
        } # fin de fila de datos / fin matriz transicion i
        diag(transition) <- 0
        mat[[i]] <- transition
      } # fin recorrer todos los datos
      
    } # fin option arc
    
    if(option == "dR-dL"){
      
      #### 
      # -------------------------
      
      mat <- list()
      n <- ncol(data)
      for (i in 1:nrow(data)){
        transition <- array(NA, dim=c(n,n))
        
        for(j in 1:n){
          if(j != n){
            for(l in (j+1):n){
              if(!is.na(data[i,j])&&!is.na(data[i,l])){
                dR <- ifelse(((data[i,j] - data[i,l])%%(2*pi))< pi, 1-cos(data[i,j] - data[i,l]), 3-cos(data[i,j] - data[i,l]-pi))
                dL <- ifelse(((data[i,j] - data[i,l])%%(2*pi))< pi, 3-cos(data[i,j] - data[i,l]-pi), 1-cos(data[i,j] - data[i,l]))
                transition[j,l] <- dL-dR
                transition[l,j] <- dR-dL       
              }
            } # recorrer cada fila
          } # if j!=n
        } # fin de fila de datos / fin matriz transicion i
        diag(transition) <- 0
        mat[[i]] <- transition
      } # fin recorrer todos los datos
      
    } # fin dR-dL
    
        
    if(option == "alpha1"){
      
      #### cos: DISTANCIAS Con cosenos   
      # -------------------------
      
      mat <- list()
      n <- ncol(data)
      for (i in 1:nrow(data)){
        transition <- array(NA, dim=c(n,n))
        
        for(j in 1:n){
          if(j != n){
            for(l in (j+1):n){
              if(!is.na(data[i,j])&&!is.na(data[i,l])){
                #coseno
                aux1 <- 1-cos(data[i,j] - data[i,l])
                transition[j,l] <- aux1
                transition[l,j] <- aux1       
              }
            } # recorrer cada fila
          } # if j!=n
        } # fin de fila de datos / fin matriz transicion i
        diag(transition) <- 0
        mat[[i]] <- transition
      } # fin recorrer todos los datos
      
    } # fin alpha1
    
    
    if(option == "alpha3"){
      
      
      # Opcion time: matriz de tiempos (segun MCUA) y distancia el coseno
      # actualizado con time2 (teniendo en cuenta los tres sectores del circulo)
      # -------------------------
      # Se usan los datos directamente (0,2pi)
      
      mat <- list()
      n <- ncol(data)
      
      for (i in 1:nrow(data)){
        times <- array(NA, dim=c(n,n))  
        
        for(j in 1:n){
          if(j != n){
            for(l in (j+1):n){
              if(!is.na(data[i,j])&&!is.na(data[i,l])){
                if((data[i,l] - data[i,j])%%(2*pi) <= (pi/2)){
                  times[j,l] <- 1-cos(data[i,j] - data[i,l])
                  times[l,j] <- 3*(1-cos(data[i,j] - data[i,l]))
                }
                if(((data[i,l] - data[i,j])%%(2*pi) > (pi/2))&&((data[i,l] - data[i,j])%%(2*pi) <= pi)){
                  times[j,l] <- 1-cos(data[i,j] - data[i,l])
                  times[l,j] <- 3-cos(data[i,j] - data[i,l]-pi)
                }
                if(((data[i,l] - data[i,j])%%(2*pi) > pi)&&((data[i,l] - data[i,j])%%(2*pi) <= (3*pi/2))){
                  times[j,l] <- 3-cos(data[i,j] - data[i,l]-pi)
                  times[l,j] <- 1-cos(data[i,j] - data[i,l])
                }
                if((data[i,l] - data[i,j])%%(2*pi) > (3*pi/2)){
                  times[j,l] <- 3*(1-cos(data[i,j] - data[i,l]))
                  times[l,j] <- 1-cos(data[i,j] - data[i,l])
                }
              }              
            } # recorrer cada fila
          } # if j!=n
        } # fin de fila de datos / fin matriz transicion i
        diag(times) <- 0
        mat[[i]] <- times
      } # fin recorrer todos los datos
      
    } # fin alpha3
    
    if(option == "alpha2"){
      
      
      # Opcion time: matriz de tiempos (segun MCUA) y distancia el coseno
      # actualizado con time2 (teniendo en cuenta los tres sectores del circulo)
      # -------------------------
      # Se usan los datos directamente (0,2pi)
      
      mat <- list()
      n <- ncol(data)
      
      for (i in 1:nrow(data)){
        times <- array(NA, dim=c(n,n))  
        
        for(j in 1:n){
          if(j != n){
            for(l in (j+1):n){
              if(!is.na(data[i,j])&&!is.na(data[i,l])){
                if((data[i,l] - data[i,j])%%(2*pi) <= (1.910633)){
                  times[j,l] <- 1-cos(data[i,j] - data[i,l])
                  times[l,j] <- 2*(1-cos(data[i,j] - data[i,l]))
                }
                if(((data[i,l] - data[i,j])%%(2*pi) > (1.910633))&&((data[i,l] - data[i,j])%%(2*pi) <= pi)){
                  times[j,l] <- 1-cos(data[i,j] - data[i,l])
                  times[l,j] <- 3-cos(data[i,j] - data[i,l]-pi)
                }
                if(((data[i,l] - data[i,j])%%(2*pi) > pi)&&((data[i,l] - data[i,j])%%(2*pi) <= (4.372552))){
                  times[j,l] <- 3-cos(data[i,j] - data[i,l]-pi)
                  times[l,j] <- 1-cos(data[i,j] - data[i,l])
                }
                if((data[i,l] - data[i,j])%%(2*pi) > (4.372552)){
                  times[j,l] <- 2*(1-cos(data[i,j] - data[i,l]))
                  times[l,j] <- 1-cos(data[i,j] - data[i,l])
                }
              }              
            } # recorrer cada fila
          } # if j!=n
        } # fin de fila de datos / fin matriz transicion i
        diag(times) <- 0
        mat[[i]] <- times
      } # fin recorrer todos los datos
      
    } # fin alpha2
    
    if(option == "alpha4"){
      
      
      # Opcion time: matriz de tiempos (segun MCUA) y distancia el coseno
      # actualizado con time2 (teniendo en cuenta los tres sectores del circulo)
      # -------------------------
      # Se usan los datos directamente (0,2pi)
      
      mat <- list()
      n <- ncol(data)
      
      for (i in 1:nrow(data)){
        times <- array(NA, dim=c(n,n))  
        
        for(j in 1:n){
          if(j != n){
            for(l in (j+1):n){
              if(!is.na(data[i,j])&&!is.na(data[i,l])){
                if((data[i,l] - data[i,j])%%(2*pi) <= (1.369438)){
                  times[j,l] <- 1-cos(data[i,j] - data[i,l])
                  times[l,j] <- 4*(1-cos(data[i,j] - data[i,l]))
                }
                if(((data[i,l] - data[i,j])%%(2*pi) > (1.369438))&&((data[i,l] - data[i,j])%%(2*pi) <= pi)){
                  times[j,l] <- 1-cos(data[i,j] - data[i,l])
                  times[l,j] <- 3-cos(data[i,j] - data[i,l]-pi)
                }
                if(((data[i,l] - data[i,j])%%(2*pi) > pi)&&((data[i,l] - data[i,j])%%(2*pi) <= (4.913747))){
                  times[j,l] <- 3-cos(data[i,j] - data[i,l]-pi)
                  times[l,j] <- 1-cos(data[i,j] - data[i,l])
                }
                if((data[i,l] - data[i,j])%%(2*pi) > (4.913747)){
                  times[j,l] <- 4*(1-cos(data[i,j] - data[i,l]))
                  times[l,j] <- 1-cos(data[i,j] - data[i,l])
                }
              }              
            } # recorrer cada fila
          } # if j!=n
        } # fin de fila de datos / fin matriz transicion i
        diag(times) <- 0
        mat[[i]] <- times
      } # fin recorrer todos los datos
      
    } # fin alpha4
    
    if(option == "alphainf"){
      
      
      # Opcion time: matriz de tiempos (segun MCUA) y distancia el coseno
      # actualizado con time2 (teniendo en cuenta los tres sectores del circulo)
      # -------------------------
      # Se usan los datos directamente (0,2pi)
      
      mat <- list()
      n <- ncol(data)
      
      for (i in 1:nrow(data)){
        times <- array(NA, dim=c(n,n))  
        
        for(j in 1:n){
          if(j != n){
            for(l in (j+1):n){
              if(!is.na(data[i,j])&&!is.na(data[i,l])){
                if(((data[i,l] - data[i,j])%%(2*pi) < pi)){
                  times[j,l] <- 1-cos(data[i,j] - data[i,l])
                  times[l,j] <- 3-cos(data[i,j] - data[i,l]-pi)
                }
                if(((data[i,l] - data[i,j])%%(2*pi) >= pi)){
                  times[j,l] <- 3-cos(data[i,j] - data[i,l]-pi)
                  times[l,j] <- 1-cos(data[i,j] - data[i,l])
                }
              }              
            } # recorrer cada fila
          } # if j!=n
        } # fin de fila de datos / fin matriz transicion i
        diag(times) <- 0
        mat[[i]] <- times
      } # fin recorrer todos los datos
      
    } # fin alpha infinito
    
    
    
    if(option == "bin"){
      # Opcion 1: matriz de ceros y unos
      # -------------------------
      # Se usan los ordenes circulares
      
      mat <- list()
      n <- ncol(cirorders)
      # para fijar el camino con ceros:
      for (i in 1:nrow(cirorders)){
        transition <- rep(1,n*n)
        transition <- matrix(transition,ncol=n)
        for(j in 1:n){
          if(j != n){transition[cirorders[i,j], cirorders[i,(j+1)]] <- 0}
          if(j==n){transition[cirorders[i,n], cirorders[i,1]] <- 0}
        }
        mat[[i]]<-transition
      }
      
    } # fin bin (antes 10)
    
    if(option == "pos"){
      
      # Opcion 2: matriz de numeros (posiciones)
      # -------------------------
      # Se usan los ordenes circulares
      
      mat <- list()
      n <- ncol(cirorders)
      # se le asigna una distancia de cero al mismo elemento
      # de 1 al siguiente en el orden segun la direccion del circulo
      # de 2 al siguiente del siguiente, y asi consecutivamente.
      for (i in 1:nrow(cirorders)){
        transition <- rep(0,n*n)
        transition <- matrix(transition,ncol=n)
        ordenation <- cirorders[i,!is.na(cirorders[i,])]
        for(j in 1:length(ordenation)){
          while(ordenation[1]!=j){
            ordenation <- c(ordenation[2:length(ordenation)],ordenation[1])
          }
          transition[j,!is.na(data[i,])] <- order(c(1:length(ordenation))[ordenation])-1
        }
        mat[[i]]<-transition
      }
      
    } # fin pos
    
    if(option == "chord"){
      
      #### 2.5: Con la cuerda
      # -------------------------
      
      mat <- list()
      n <- ncol(data)
      for (i in 1:nrow(data)){
        transition <- array(NA, dim=c(n,n))
        
        for(j in 1:n){
          if(j != n){
            for(l in (j+1):n){
              if(!is.na(data[i,j])&&!is.na(data[i,l])){
                # calculando la cuerda, por el teorema del coseno
                aux1 <- sqrt(2 - 2*cos(data[i,j] - data[i,l])) 
                # calculando la cuerda, por el seno
                #aux1 <- 2*abs(sin((data[i,j]-data[i,l])/2))
                transition[j,l] <- aux1
                transition[l,j] <- aux1
              }
            } # recorrer cada fila
          } # if j!=n
        } # fin de fila de datos / fin matriz transicion i
        diag(transition) <- 0
        mat[[i]] <- transition
      } # fin recorrer todos los datos
      
    } # fin chord
    
    
    if(option == "time"){
      
      
      # Opcion 3.1: matriz de tiempos (segun MCUA) y distancia el coseno
      # -------------------------
      # Se usan los datos directamente (0,2pi)
      
      mat <- list()
      n <- ncol(data)
      
      for (i in 1:nrow(data)){
        rotation <- TRUE
        times <- array(NA, dim=c(n,n))  
        
        for(j in 1:n){
          if(j != n){
            for(l in (j+1):n){
              if(!is.na(data[i,j])&&!is.na(data[i,l])){
                
                # el coseno entre cada par de dos puntos:
                aux1 <- 1-cos(data[i,j] - data[i,l])
                
                if(data[i,j] <= pi){
                  rotation <- ifelse((data[i,l] - data[i,j])%%(2*pi) <= pi, TRUE, FALSE)
                }
                if(data[i,j] > pi){
                  rotation <- ifelse((data[i,j] - data[i,l])%%(2*pi) <= pi, FALSE, TRUE)
                }
                
                if(rotation){
                  times[j,l] <- aux1*(2/(3*pi))
                  times[l,j] <- aux1*(2/pi)          
                }
                if(!rotation){
                  times[j,l]<-aux1*(2/pi)
                  times[l,j]<-aux1*(2/(3*pi))
                }
                
              }
              
            } # recorrer cada fila
          } # if j!=n
        } # fin de fila de datos / fin matriz transicion i
        diag(times) <- 0
        mat[[i]] <- times
      } # fin recorrer todos los datos
      
    } # fin time with cosine
    

    
    
    
    ###################
    # Unimos en una sola matriz a las matrices que representan cada orden 

      n <- ncol(mat[[1]])
      matTSP <- matrix(ncol=n, nrow=n)
      if(length(mat)==1){matTSP <- ws[1]*mat[[1]]}
      if(length(mat)>1){
        aux <- NULL
        for(i in 1:n){
          for(j in 1:n){
            for(l in 1:length(mat)){
              aux[l] <- mat[[l]][i,j]
            }
            matTSP[i,j] <- weighted.mean(aux, ws, na.rm = TRUE)
          }
        }
      }
    
    return(list(mat = mat, matTSP = matTSP))
    
  } # fin CORAM
  
  hodgefusion <- function(data, ws=NULL, control.method=c("bin", "pos", "cos","cos-", "cmean", "mrl", "e3", "dR-dL", "ave","qua","nat","natb")){
    
    tita<-data
    if(!is.matrix(tita)&is.vector(tita)){tita<-rbind(tita)}
    control.method <- match.arg(control.method)
    if(is.null(ws)){ws<-rep(1,nrow(tita))}
    n <- ncol(tita)
    X <- array(0,dim=c(n,n))
    
    
    if(control.method!="dR-dL" && control.method!="nat"&& control.method!="natb"){
    Phi <- list(array(0,dim=c(n,n,n)))
    length(Phi) <- nrow(tita)
    if(nrow(tita)>1){
      for(i in 1:(nrow(tita)-1)){Phi[[i+1]] <- array(0,dim=c(n,n,n))}
    }  
    Phimedia <- array(0,dim=c(n,n,n))
    
    if(control.method=="i0"){
      data_z <- array(dim=dim(tita))
      for (i in 1:nrow(tita)){data_z[i,] <- as.numeric((tita[i,]-tita[i,1]))%%(2*pi)}
      tita <- data_z
    }
    
    for(a in 1:nrow(tita)){
      Phi[[a]] <- array(0,dim=c(n,n,n))
      if(control.method=="pos"){
        # posiciones de fuera menos posiciones de dentro
        P<-diag(x = 0, n, n)
        for(i in 1:n){
          for(j in i:n){
            if(i!=j){
              if(all(!is.na(c(tita[a,i], tita[a,j])))){
              if(tita[a,i]-tita[a,j]>0){
                Fij <- sum((tita[a,]>tita[a,j])&(tita[a,]<tita[a,i]))
                Dij <- sum((tita[a,]>tita[a,i])|(tita[a,]<tita[a,j]))
                P[i,j] <- Fij-Dij
              }
              if(tita[a,i]-tita[a,j]<0){
                Fij <- sum((tita[a,]>tita[a,j])|(tita[a,]<tita[a,i]))
                Dij <- sum((tita[a,]>tita[a,i])&(tita[a,]<tita[a,j]))
                P[i,j] <- Fij-Dij               
              }
              if(tita[a,i]-tita[a,j]==0){P[i,j]<-0}              
            }
            }
            if(any(is.na(c(tita[a,i], tita[a,j])))){P[i,j]<-NA}
          }
        }
        #P[lower.tri(P, diag=FALSE)]<- -P[upper.tri(P, diag=FALSE)] no los coloca bien
        P <- P-t(P)
      } 
      
      if(control.method=="e3"){E<-CORAM(rbind(tita[a,]), option = "alpha3", ws = ws)$mat[[1]]}
      for(i in 1:n){
        for(j in 1:n){
          for(k in 1:n){
            sigijk <- sign(tita[a,j]-tita[a,i])+sign(tita[a,k]-tita[a,j])+sign(tita[a,i]-tita[a,k])
            if(control.method=="bin"){
              # solo los signos
              Phi[[a]][i,j,k] <- sigijk
            }
            if(control.method=="pos"){              
              Phi[[a]][i,j,k] <- sigijk*abs(P[i,j]+P[j,k]+P[k,i])
            }
            if(control.method=="cos"){    
              # manteniendo el sign del triplete y cos entre ijk
              Phi[[a]][i,j,k] <- sigijk*((1+cos(tita[a,j]-tita[a,i]))+
                                           (1+cos(tita[a,k]-tita[a,j]))+(1+cos(tita[a,i]-tita[a,k])))
            }
            if(control.method=="qua"){   
              # escribir a, b y c:
              aux <- c(tita[a,i],tita[a,j],tita[a,k])
              b <- (aux[2]-aux[1])
              c <- (aux[3]-aux[2])
              d <- (2*pi)-(b+c)
              Phi[[a]][i,j,k] <- sigijk*(1-(((b/(2*pi))^2)+((c/(2*pi))^2)+((d/(2*pi))^2)))
            }
#             if(control.method=="cos-"){    
#               # manteniendo el sign del triplete y cos entre ijk
#               Phi[[a]][i,j,k] <- sigijk*((1-cos(tita[a,j]-tita[a,i]))+
#                                            (1-cos(tita[a,k]-tita[a,j]))+(1-cos(tita[a,i]-tita[a,k])))
#             }
            if(control.method=="cmean"){
              # cos entre cada uno y la media de ellos
              cmeanijk <- cirmean(c(tita[a,i],tita[a,j],tita[a,k]))
              Phi[[a]][i,j,k] <- sigijk*((cos(tita[a,j]-cmeanijk))+
                                           (cos(tita[a,k]-cmeanijk))+(cos(tita[a,i]-cmeanijk)))
            }
            if(control.method=="ave"){
              aux1 <- min((tita[a,j] - tita[a,i])%%(2*pi), (tita[a,i] - tita[a,j])%%(2*pi))
              aux2 <- min((tita[a,k] - tita[a,j])%%(2*pi), (tita[a,j] - tita[a,k])%%(2*pi))
              aux3 <- min((tita[a,i] - tita[a,k])%%(2*pi), (tita[a,k] - tita[a,i])%%(2*pi))
              Phi[[a]][i,j,k] <- sigijk*cirmean(c(aux1, aux2, aux3))
              
            }
#           if(control.method=="smean"){
#               # sen entre cada uno y la media de ellos ### NO DEFINE CORRECTAMENTE
#               cmeanijk <- cirmean(c(tita[a,i],tita[a,j],tita[a,k]))
#               Phi[[a]][i,j,k] <- sigijk*((sin(tita[a,j]-cmeanijk))+
#                                            (sin(tita[a,k]-cmeanijk))+(sin(tita[a,i]-cmeanijk)))
#           }
            if(control.method=="mrl"){
              Phi[[a]][i,j,k] <- sigijk*(1-mrl(t(c(tita[a,i],tita[a,j],tita[a,k]))))
              
            }
#           if(control.method=="marc"){
#               Phi[[a]][i,j,k] <- sigijk*min(c((tita[a,j]-tita[a,i])%%(2*pi),(tita[a,k]-tita[a,j])%%(2*pi),(tita[a,i]-tita[a,k])%%(2*pi)))            
#           }
            if(control.method=="e3"){
              Phi[[a]][i,j,k] <- sigijk*(E[i,j]+E[j,i]+E[j,k]+E[k,j]+E[k,i]+E[i,k])
            }
            if(control.method=="i0" && i==1 &&j!=1 && k!=1){
              Phi[[a]][i,j,k] <- tita[a,k] - tita[a,j]  
              Phi[[a]][k,1,j] <- tita[a,k] - tita[a,j]  
              Phi[[a]][j,k,1] <- tita[a,k] - tita[a,j]  
            }
          }      
        }     
      }  
    }
   
    aux <- rep(0,nrow(tita))
    for(i in 1:n){
      for(j in 1:n){
        for(k in 1:n){
          for(a in 1:nrow(tita)){
            aux[a] <- Phi[[a]][i,j,k]
          }
          if(control.method!="i0"){Phimedia[i,j,k] <- weighted.mean(aux, ws, na.rm=TRUE)}
          if(control.method=="i0"){Phimedia[i,j,k] <- cirmean(aux, ws)}
          X[i,j] <- X[i,j]+Phimedia[i,j,k]
        }
      }
    }
    
    }
    
   
    if(control.method=="dR-dL"){X <- CORAM(rbind(data), option = control.method, ws = ws)$matTSP}
    #quitando uno de los que tienen valor max:
    #out <- ceiling(which.max(X)/nrow(X))


    if((control.method!="nat")&&(control.method!="natb")){
    # quitando el que suma maximo en valor absoluto:
    out<-which.max(apply((X^2),1,sum))
    X1 <- X[-out,-out]
    colnames(X1) <- c(1:n)[-out]
    rownames(X1) <- c(1:n)[-out]
    scores <- -apply(X1, 1, mean)
    dxly <- (sum(X1^2))+(2*sum(X[out,-out]^2)/(n-1)) # distancia(X[-l],Y)
    tita_order <-  c(as.numeric(names(sort(scores))),out)
    #sce<-msce(tita, order(tita_order), ws)$msce
#     tita_order
#     sce
#     data<-cbind(data, data[,1])
#     data<-data[,-1]
#     tita<-data
    #tau <- mcirktau(tita, order(tita_order), ws)$mtau
    Xe <- matrix(ncol=(n-1), nrow=(n-1))
    sest <- -X[out,-out]/(n-1)
    names(sest) <- names(scores)
    for(i in 1:(n-1)){
      for(j in i:(n-1)){
        Xe[i,j] <- sest[j]-sest [i]
        Xe[j,i] <- sest[i]-sest [j]
      }
    }
    scoresest <-  -apply(Xe,1,mean)
    names(scoresest) <- names(scores)
    c(as.numeric(names(sort(scoresest))),out) # comprobando que tiene el mismo orden Xe (X[-l,-l] en el papel) que X1 (Y[-l,-l] en el papel).
    error1<-norm(X1-Xe, type="2")/norm(X1,type="2")
    error2<-norm(Xe, type="2")/norm(X1,type="2")
    #resultado <-  list(aggre_order=tita_order, msce=sce, mtau=tau, scores=scores, out=out, error1=error1, error2=error2)
    resultado <-  list(aggre_order=tita_order, scores=scores, out=out, error1=error1, error2=error2)
    }

if(control.method=="natb"){
  
  resultados <- list()
  for(i0 in 1:n){
    Phimedia <- array(0,dim=c(n,n,n)) # agregada
    X <- array(0,dim=c(n,n)) # info agregada
    Phi <- list(array(0,dim=c(n,n,n)))
    length(Phi) <- nrow(data)
    rot <-matrix(ncol=n,nrow=nrow(data))
    mat <- list()
    for(a in 1:nrow(rot)){
      Phi[[a]] <- array(0,dim=c(n,n,n))
      rot[a,] <- (data[a,]-data[a,i0])%%(2*pi)
      posi <- order(rot[a,!is.na(rot[a,])])
      transition <- rep(0,n*n)
      transition <- matrix(transition,ncol=n)
      for(j in 1:n){
        if(j != n){
          transition[posi[j], posi[(j+1)]] <- 1
          transition[posi[j], posi[(j-1)]] <- 1
        }
        if(j==n){
          transition[posi[n], posi[1]] <- 1
          transition[posi[j], posi[(j-1)]] <- 1
        }
        if(j==1){transition[posi[1], posi[n]] <- 1}
      }
      mat[[a]]<-transition 
      for(i in 1:n){
        for(h in 1:n){
          for(k in 1:n){
            sigijk <- sign(rot[a,j]-rot[a,i])+sign(rot[a,k]-rot[a,j])+sign(rot[a,i]-rot[a,k])
            if(i==i0&&h!=i0&&k!=i0){Phi[[a]][i,h,k] <- sigijk*(mat[[a]][h,k])}
            if(h==i0&&i!=i0&&k!=i0){Phi[[a]][i,h,k] <- sigijk*(mat[[a]][k,i])}
            if(k==i0&&h!=i0&&i!=i0){Phi[[a]][i,h,k] <- sigijk*(mat[[a]][i,h])}
          }
        }
      }
      
    }  
    
    aux <- rep(0,nrow(rot))
    for(i in 1:n){
      for(j in 1:n){
        for(k in 1:n){
          for(a in 1:nrow(rot)){
            aux[a] <- Phi[[a]][i,j,k]
          }
          Phimedia[i,j,k] <- weighted.mean(aux, ws, na.rm=TRUE)
          X[i,j] <- X[i,j]+Phimedia[i,j,k]
        }
      }
    }
    
#     matTSP <- matrix(ncol=n, nrow=n)
#     if(length(mat)==1){matTSP <- ws[1]*mat[[1]]}
#     if(length(mat)>1){
#       aux <- NULL
#       for(i in 1:n){
#         for(j in 1:n){
#           for(l in 1:length(mat)){
#             aux[l] <- mat[[l]][i,j]
#           }
#           matTSP[i,j] <- weighted.mean(aux, ws, na.rm = TRUE)
#         }
#       }
#     }
#     for(i in 1:n){
#       for(h in 1:n){
#         for(k in 1:n){
#           if(i==i0&&h!=i0&&k!=i0){Phimedia[i,h,k] <- n*(matTSP[h,k])}
#           if(h==i0&&i!=i0&&k!=i0){Phimedia[i,h,k] <- n*(matTSP[k,i])}
#           if(k==i0&&h!=i0&&i!=i0){Phimedia[i,h,k] <- n*(matTSP[i,h])}
#           X[i,h] <- X[i,h]+Phimedia[i,h,k]
#         }
#       }
#     }
     resultados[[i0]] <- X
  } # fin de para cada i0   
normas <- NULL
for(i in 1:length(resultados)){
  normas[i] <- sum(X^2)
}
X <- resultados[[which.max(normas[i])]]
out<-which.max(apply((X^2),1,sum))
X1 <- X[-out,-out]
colnames(X1) <- c(1:n)[-out]
rownames(X1) <- c(1:n)[-out]
scores <- -apply(X1, 1, mean)
tita_order <-  c(as.numeric(names(sort(scores))),out)
resultado <-  list(aggre_order=tita_order)
} # fin nat binary
    
    if(control.method=="nat"){
      resultados <- list()
    for(i0 in 1:n){
    Phimedia <- array(0,dim=c(n,n,n)) # agregada
    X <- array(0,dim=c(n,n)) # info agregada
    rot <-matrix(ncol=n,nrow=nrow(data))
    for(j in 1:nrow(data)){ # recorre experimentos
        rot[j,] <- (data[j,]-data[j,i0])%%(2*pi)
    }
    avecir <- apply(rot, 2, cirmean, ws)
        for(i in 1:n){
          for(h in 1:n){
            for(k in 1:n){
              if(i==i0&&h!=i0&&k!=i0){Phimedia[i,h,k] <- n*(avecir[k]-avecir[h])}
              if(h==i0&&i!=i0&&k!=i0){Phimedia[i,h,k] <- n*(avecir[i]-avecir[k])}
              if(k==i0&&h!=i0&&i!=i0){Phimedia[i,h,k] <- n*(avecir[h]-avecir[i])}
              X[i,h] <- X[i,h]+Phimedia[i,h,k]
            }
          }
        }
        # se colova Phimedia^j 
        #Phimedia[[i0]] <- Phimedia
    resultados[[i0]] <- X
    } # fin de para cada i0   
    normas <- NULL
    for(i in 1:length(resultados)){
      normas[i] <- sum(X^2)
    }
    X <- resultados[[which.max(normas[i])]]
    out<-which.max(apply((X^2),1,sum))
    X1 <- X[-out,-out]
    colnames(X1) <- c(1:n)[-out]
    rownames(X1) <- c(1:n)[-out]
    scores <- -apply(X1, 1, mean)
    tita_order <-  c(as.numeric(names(sort(scores))),out)
    resultado <-  list(aggre_order=tita_order)
    } # fin de CHnat
    return(resultado)  
  } # fin Hodge Fusion
  
  cirmean<-function (data, ws=NULL) {
    # data is a vector
    if(is.null(ws)){ws<-rep(1,length(data))}
    if(any(is.na(data))){
      ws<-ws[complete.cases(data)]
      ws <- ws/sum(ws)}
    data <- data[complete.cases(data)]
    A <- sum(ws*sin(data))
    B <- sum(ws*cos(data))
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
    if((A==0)&&(B==0)){result <- 0}
    return(result)
  }
  
  rotation <- function(ordenation){
    # rota el orden introducido en el argumento "ordenation" 
    # hasta que el primer elemento sea el 1
    while(ordenation[1] != 1){
      ordenation <- c(ordenation[2:length(ordenation)], ordenation[1])
    }
    return(ordenation)
  }
  
  bordacircular <- function(data, ws, control.method=c("pos", "cirmean", "cirmed")){
    
  if(control.method=="pos"){
    ###################### Posiciones ######################    
    MC_order <- matrix(nrow=ncol(data), ncol=ncol(data))
    for(j in 1:ncol(data)){
      data_z <- array(dim=dim(data))
      for (i in 1:nrow(data)){data_z[i,] <- as.numeric((data[i,]-data[i,j]))%%(2*pi)}
      cirorders <- t(apply(data_z, 1, order))
      cirorders <- t(apply(cirorders, 1, rotation))
      posorders <- t(apply(cirorders, 1, order))
      meanord <- apply((2*pi/ncol(data))*(posorders-1),2,cirmean,ws)
      MC_order[j,] <- order(meanord)
    }
    MC_order <- t(apply(MC_order, 1, rotation))
    SOL1 <- rbind(unique(rbind(MC_order), MARGIN=1))   
    msces1 <- array(dim=nrow(SOL1))
    for(i in 1:nrow(SOL1)){msces1[i]<-msce(data, order(SOL1[i,]), ws=ws)$msce}
    aggre_order <- SOL1[which.min(msces1),]
    sce <- min(msces1)
    tau <- mcirktau(data, order(aggre_order), ws)$mtau
  }
    
    if(control.method=="cirmean" || control.method=="cirmed"){
      
    ###################### Medias o Medianas Circulares ######################  
    MC_order <- matrix(nrow=ncol(data), ncol=ncol(data))
    for(j in 1:ncol(data)){
      data_z <- array(dim=dim(data))
      for (i in 1:nrow(data)){data_z[i,] <- as.numeric((data[i,]-data[i,j]))%%(2*pi)}
      if(control.method=="cirmean"){MC_0 <- apply(data_z, 2, cirmean, ws)}
      if(control.method=="cirmed"){
        MC_0 <- rep(NA, ncol(data))
        for(k in 1:ncol(data)){
          MC_0[k] <- suppressWarnings(medianCircular(data_z[,k], na.rm=TRUE))
        }
      }
      MC_order[j,] <- order(MC_0)
    }
    
    MC_order <- t(apply(MC_order, 1, rotation))
   SOL1 <- rbind(unique(rbind(MC_order), MARGIN=1))   
      msces1 <- array(dim=nrow(SOL1))
      for(i in 1:nrow(SOL1)){
        msces1[i]<-msce(data, order(SOL1[i,]), ws=ws)$msce
      }
      
      aggre_order <- SOL1[which.min(msces1),]
      sce <- min(msces1)
    tau <- mcirktau(data, order(aggre_order), ws)$mtau
    
    }
  
  return(list(aggre_order=aggre_order, msce=sce, mtau=tau))
      
        
  }
  
  TSPmethod <- function (data, ws, coef, control.method){
    
    
    distmatT <- CORAM(data, option = control.method, ws = ws)
    matTSPT<-as.ATSP(distmatT$matTSP)
    pmethods <- c("farthest_insertion", "nearest_insertion", "nn",
                  "cheapest_insertion", "arbitrary_insertion", "repetitive_nn")
    veces <- ifelse(coef<1,1,coef)
    methods <- c(rep(pmethods, veces*ncol(data)))
    # los repitimos porque son muy inestables y susceptibles de sacar diferentes 
    # resultados en cada ejecucion
    toursT <- list()
    length(toursT) <- length(methods)
    for(i in 1:length(methods)){
      toursT[[i]] <- solve_TSP(matTSPT, method = methods[i])
    }
    solutionT <- matrix(ncol=(ncol(data)+1), nrow=length(toursT))
    for(i in 1:length(toursT)){
      solutionT[i, c(1:ncol(data))] <- as.numeric(labels(toursT[[i]]))
      solutionT[i, ncol(data)+1] <- attr(toursT[[i]],"tour_length")    
    }
  
#     
#     
#     if(length(control.method)==2){
#       distmatT <- CORAM(data, option = control.method[2], ws = ws)
#       matTSPT<-as.ATSP(distmatT$matTSP)
#       toursT <- list()
#       length(toursT) <- length(methods)
#       for(i in 1:length(methods)){
#         toursT[[i]] <- solve_TSP(matTSPT, method = methods[i])
#       }
#       solutionT2 <- matrix(ncol=(ncol(data)+1), nrow=length(toursT))
#       for(i in 1:length(toursT)){
#         solutionT2[i, c(1:ncol(data))] <- as.numeric(labels(toursT[[i]]))
#         solutionT2[i, ncol(data)+1] <- attr(toursT[[i]],"tour_length")    
#       }
#       solutionT <- rbind(solutionT, solutionT2)
#     } # fin si 2 control.method   
#     
    
    ## en solution tenemos todas las soluciones para los metodos de arriba
    if(length(control.method)==1){rownames(solutionT) <- methods}
   # if(length(control.method)==2){rownames(solutionT) <- rep(methods,2)}
    solution <- rbind(unique(rbind(solutionT), margin=1))
    
    # nos quedamos con n los mejores en cuanto a longitud de tour (si no hay mas de n, con los que haya)
    fin <- floor(coef*ncol(data))
    if(fin==0){fin <- 1}
    if(nrow(solution)<fin){fin <- nrow(solution)}
    solution3 <- matrix(ncol=ncol(data),nrow=fin) 
    for(i in 1:fin){
      index <- which(round(as.numeric(solution[,(ncol(data)+1)]),6)==round(as.numeric(levels(as.factor(solution[,(ncol(data)+1)])))[i],6))[1]
      solution3[i,] <- solution[index,c(1:ncol(data))]   
    }
    
    solution3 <- rbind(solution3[!is.na(solution3[,1]),])
#     if(control.method=="alpha1" || control.method=="chord"){
#       solution3 <- rbind(solution3, solution3[1,c(ncol(data):1)])
#     }
    solution3 <- t(apply(rbind(solution3), 1, rotation))
    solution3 <- rbind(unique(rbind(solution3), MARGIN=1))
  
  SOL1 <- solution3
if(coef>0){
  msces1 <- array(dim=nrow(SOL1))
  for(i in 1:nrow(SOL1)){
    msces1[i]<-msce(data, order(SOL1[i,]), ws=ws)$msce
  }
  
  aggre_order <- SOL1[which.min(msces1),]
  sce <- min(msces1)
  tau <- mcirktau(data, order(aggre_order), ws)$mtau
  
  minlength <- min(solutionT[,(ncol(data)+1)])
  elmin <- which(solutionT[,(ncol(data)+1)]==minlength)[1]
  mtour <- solutionT[elmin, 1:ncol(data)]
  mt_msce <- msce(data, order(mtour), ws)$msce
  X <- distmatT$matTSP
  scores <- array(,dim=c(length(aggre_order)))
  for(i in 1:length(aggre_order)){
    if(i<length(aggre_order)){sig<-aggre_order[i+1]}
    if(i==length(aggre_order)){sig <- aggre_order[1]}
    scores[i] <- X[aggre_order[i], sig]
  }
  resultado <- list(aggre_order=aggre_order, msce=sce, mtau=tau, mintour=mtour, mt_msce=mt_msce, tour_length=minlength, scores=scores)
  }
if(coef==0){
  aggre_order <- SOL1[1,]
  resultado <- list(aggre_order=aggre_order)
}
  return(resultado)
} 
  
  
  
  ###################### control of arguments ######################
  
  if(is.data.frame(data)){data <- as.matrix(data)}
  if(!(is.matrix(data))){stop("The argument data must be in matrix format")}
  if(!nrow(data)>1){stop("There is no possible fusion with one row in data")}
  if(!ncol(data>2)){stop("At least 3 elements are needed to do circular fusion")}
  if(any(apply(is.na(data),1,all))){data <- data[-which(apply(is.na(data),1,all)),]}
  if(any(apply(is.na(data),2,all))){data <- data[,-which(apply(is.na(data),2,all))]}
  if(is.null(ws)){ws<-(rep(1,nrow(data)))/nrow(data)}
  if(sum(ws)!=1){ws<-ws/sum(ws)}
  method <- match.arg(method, c("Naive", "CB", "CMC", "TSP", "CH"))
  
  if(method=="Naive"){control.method <- match.arg(control.method, c("tau", "MSCE"))}
  if(method=="CB"){control.method <- match.arg(control.method, c("pos", "cirmean", "cirmed"))}
  if(method=="CMC"){control.method <- match.arg(control.method, c("1", "2", "3", "4m", "4c"))}
  if(method=="TSP"){
    control.method[1] <- match.arg(control.method[1], c("bin", "pos", "alpha1", "alpha2", "alpha3", "alpha4", "alphainf", "time", "arc", "chord"))
#     if(length(control.method)==2){
#       control.method[2] <- match.arg(control.method[2], c("bin", "pos", "alpha1", "alpha2", "alpha3", "alpha4", "alphainf", "time", "arc", "chord"))
#     }
  }
  if(method=="CH"){control.method <- match.arg(control.method, c("bin", "pos", "cos", "cmean", "mrl", "e3", "dR-dL", "ave", "qua","nat", "natp", "natb"))}
  
  # control of path between all points
  if(nrow(which(is.na(data),arr.ind = TRUE))>1){
    a<-cbind(combn(as.numeric(names(table(which(is.na(data[,]),arr.ind = TRUE)[,2]))),2))
    for(j in 1:ncol(a)){
      if(length(which(!is.na(apply(data[,a[,j]],1,sum))))==0){
        stop(paste("There is no path between points",a[1,j],"and",a[2,j] ))
      }
    }    
  }

if(method=="Naive"){
  output <- circularOrderFussion(data, ws, type=0, control.checking=control.method)
}

if(method=="CB"){
  output <- bordacircular(data, ws, control.method=control.method)
}

if(method=="CMC"){
  if(control.method!="4m" && control.method!="4c"){
    output <- circularOrderFussion(data, ws, type=as.numeric(control.method))
  }
  if(control.method=="4m"){
    output <- circularOrderFussion(data, ws, type=4, control.checking="majority")
  }
  if(control.method=="4c"){
    output <- circularOrderFussion(data, ws, type=4, control.checking="count")
  }
}

if(method=="TSP"){
  output <- TSPmethod(data, ws, control.method=control.method, coef)  
}
  
  ###################### Circular Hodge Theory ######################
  
  if(method=="CH"){
    if(control.method=="natp"){
      positions <- matrix(ncol=ncol(data),nrow=nrow(data))
      for(i in 1:nrow(data)){
         posi <- order(order(data[i,!is.na(data[i,])]))
         positions[i,!is.na(data[i,])] <- posi*(2*pi/length(posi))
      }
      control.method <- "nat"
      output <- hodgefusion(positions, control.method=control.method,  ws=ws)
    }
    else{output <- hodgefusion(data, control.method=control.method,  ws=ws)}
  }  
  
  return(output)
  
} # end ACO
