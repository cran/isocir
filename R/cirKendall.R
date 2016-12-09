"cirKendall" <- function(phi1, phi2, test = FALSE, control.test = c("noteq", "upper", "lower")){
  
  # DATE WRITTEN: 29 Sep 2011          LAST REVISED:  27 Ene 2013
  # AUTHOR: Sandra Barragan.
  # DESCRIPTION: Computes the Circular Kendall Tau.
  # REFERENCE: Fisher, N.I. (1993) Statistical analysis of circular data.
  
  # control de los argumentos:
  if (NROW(phi1)!=NROW(phi2)){stop("phi1 and phi2 must have the same number of observations")}
  control.test <- match.arg(control.test) 
  if(any(is.na(phi1))||any(is.na(phi2))){
    nNA1 <- ifelse(complete.cases(phi1), 0, 1)
    if(sum(nNA1) >= 1){
      arg1 <- phi1[nNA1 == 0]
      arg2 <- phi2[nNA1 == 0]
      nNA2 <- ifelse(complete.cases(arg2), 0, 1)
      if(sum(nNA2) >= 1){
        phi2 <- arg2[nNA2 == 0]
        phi1 <- arg1[nNA2 == 0]
      }
      if(sum(nNA2)==0){
        phi1 <- arg1
        phi2 <- arg2
      }       
    }
    if(sum(nNA1)==0){
      phi1 <- phi1[complete.cases(phi2)]
      phi2 <- phi2[complete.cases(phi2)]
    }
  }
  # control.test tiene la hipotesis alternativa del contraste a realizar si test=TRUE
  
  n <- length(phi1)  
  pairs <- cbind(phi1, phi2)  
  combs <- as.matrix(combn(1:length(phi1),3))
  aux1 <- matrix()
  delta <- vector()
  # se calculan los tripletes concordantes y discordantes
  for (i in 1:ncol(combs)){
    aux1 <- pairs[combs[,i],]
    delta[i] <- prod(sign(aux1[1,] - aux1[2,])*sign(aux1[2,] - aux1[3,])*sign(aux1[3,] - aux1[1,]))
  }
  DELTA <- sum(delta/(length(delta) - length(delta[delta == 0])))
  y <- n*DELTA
  result <- DELTA
  #  hasta aqui sin test=FALSE y la funcion devolvera un unico valor: la Tau Circular de Kendall entre phi1 y phi2
  
  if (test == TRUE){
    if (n < 8){
      # tabulation of the null hypothesis for n<8
      matrix3 <- rbind(c(3,-3),c(0.5,0.5))
      matrix4 <- rbind(c(4,0,-4),c(1/6,4/6,1/6))
      matrix5 <- rbind(c(5,2,1,0,-1,-2,-5),c(1/24,5/24,5/24,2/24,5/24,5/24,1/24))
      matrix6 <- rbind(c(6,3.6,2.4,1.2,0,-1.2,-2.4,-3.6,-6),c(1/120,6/120,12/120,23/120,36/120,23/120,12/120,6/120,1/120))
      matrix7 <- rbind(c(7,5,3.8,3.4,2.6,2.2,1.8,1.4,1,0.6,0.2,-0.2,-0.6,-1,-1.4,-1.8,-2.2,-2.6,-3.4,-3.8,-5,-7.5),(1/720)*c(1,7,14,21,14,21,63,44,28,70,77,77,70,28,44,63,21,14,21,14,7,1))
      table12 <- list(matrix3,matrix4,matrix5,matrix6,matrix7)
      auxp <- table12[[(n-2)]]	
      if (control.test == "noteq"){
        pvalue <- 2*sum(auxp[2, c(1:ncol(auxp))[auxp[1,] > abs(y)]])
      }
      if (control.test == "upper"){
        pvalue <- sum(auxp[2, c(1:ncol(auxp))[auxp[1,] > y]])
      }
      if (control.test == "lower"){
        pvalue <- sum(auxp[2, c(1:ncol(auxp))[auxp[1,] < y]])
      }
    } # end n<8
    if (n >= 8){
      message <- NULL
      alfa1 <- c(0.1, 0.05, 0.025, 0.01, 0.005, 0.001)
      table12b <- rbind(c(8, 2.1, 2.78, 3.4, 4.25, 4.74, 6.23),c(9, 2.07, 2.72, 3.33, 4.16, 4.64, 6.08), c(10, 2.04, 2.68, 3.28, 4.09, 4.56, 5.96),c(11, 2.01, 2.65, 3.24, 4.02, 4.5, 5.85), c(12, 1.99, 2.62, 3.2, 3.97, 4.44, 5.77), c(13, 1.97, 2.59, 3.17, 3.93, 4.4, 5.7), c(14, 1.96, 2.57, 3.14, 3.89, 4.36, 5.64), c(15, 1.95, 2.56, 3.12, 3.86, 4.33, 5.59))
      if(n <= 15){
        values <- table12b[c(1:nrow(table12b))[table12b[,1]==n],2:ncol(table12b)]
        if (y < values[1]){pvalue <- 0.1
                           message <- "pvalue more than 0.1"}
        if (y > values[6]){pvalue <- 0.001
                           message <- "pvalue lower than 0.001"}
        if (y >= values[1] && y <= values[6]){
          pvalue <- c(alfa1[max(c(1:length(values))[values <= y])],alfa1[min(c(1:length(values))[values >= y])])
        }
      } # end n<=15
      if (n > 15){
        values <- c(1.77, 2.31, 2.81, 3.42, 3.85, 4.85)
        if (y < values[1]){pvalue <- 0.1
                           message <- "pvalue more than 0.1"}
        if (y > values[6]){pvalue <- 0.001
                           message <- "pvalue lower than 0.001"}
        if (y >= values[1] && y <= values[6]){
          pvalue <- c(alfa1[max(c(1:length(values))[values <= y])],alfa1[min(c(1:length(values))[values >= y])])
        }
      } # end n>15
      if (control.test == "noteq"){
        message <- ifelse(pvalue==0.1,"pvalue more than 0.2",message)
        message <- ifelse(pvalue==0.001,"pvalue lower than 0.002",message)
        pvalue <- 2*pvalue
      }
      if (control.test == "upper"){
        pvalue <- pvalue
      }
      if (control.test == "lower"){
        pvalue <- pvalue
      }
    } # end n>=8
    result <- list(DELTA = result, pvalue = pvalue, message)
    # el resultado de la funcion si test=TRUE
    # es una lista con la Tau Circular de Kendall en el elemento $DELTA
    # el pvalor del test en $pvalue y un mensaje segun pvalue
    
  } # end test
  
  return(result)
} # fin cirKendall
