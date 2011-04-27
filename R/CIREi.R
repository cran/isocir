"CIREi"<-function(data, levels=c(1:nrow(data)), isotropic=TRUE, graphic=FALSE, stack=TRUE,...){

# DATE WRITTEN: 04 Mar 2010          LAST REVISED:  04 Abr 2010
# AUTHOR: Sandra Barragan based on the SAS routines written by Miguel A. Fernandez.
# DESCRIPTION: Computes the algorithm to obtain the CIRE.
# REFERENCE: Rueda, C., Fernandez, M. A. and Shyamal, P. (2009). 
#             Estimation of parameters subject to order restrictions 
#             on a circle with application to estimation of phase angles of cell-cycle genes. 
#             JASA.
# SEE ALSO: cirPAVA, cirmean, cirSCE.


mmeans<-function(data){
means<- matrix(0,nrow=length(data),ncol=length(data))
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

cirPAVA<-function(data){
means<-mmeans(data)
data1<-data
problems<-ifelse(data[1:length(data)-1]-data[2:length(data)]>0,1,0)
while(sum(problems)!=0){
 pos<-list()
 k<-1
 for(i in 1:length(problems)){
  if(problems[i]==1){
   if((i!=1)&&problems[i-1]==0){pos[[k]]<-vector()}
    if(i==1){pos[[k]]<-vector()}
    pos[[k]]<-c(pos[[k]],c(1:length(problems))[i])
    }
  if((problems[i]==0) && (i!=1)){
   if((problems[i-1]==1)&&(i!=length(problems))){
    k<-k+1
    }
   }
  } # close for(i)
 for(j in 1:length(pos)){
  a<-min(pos[[j]])
  a<-min(c(1:length(data))[data==data[a]])
  b<-max(pos[[j]])+1
  b<-max(c(1:length(data))[data==data[b]])
  data[a:b]<-rep(means[a,b],(b-a)+1)
  } # close for(j)
 problems<-ifelse(data[1:length(data)-1]-data[2:length(data)]>0,1,0)
 } # close while
SCE<-sum(1-cos(data1-data))
return(list(data=data,SCE=SCE))
}

funCIRE<-function(data,means,mrl){

data1<-data
q<-length(data)
lA1<-list()
lA2<-list()
lA3<-list()

### G1
for(s in 0:q){
 G1<-array()
 A1<-array()
 if(s!=0){
  G1<-data[1:s]
  A1<-rbind(G1)
  aux<-0
  for(i in s:1){
   fin<-nrow(A1)
   for(j in 1:fin){
    if(A1[j,i]<pi){
     solution<-A1[j,]
     A1<-rbind(A1,solution)
     }
    if(A1[j,i]<(2*pi)&&A1[j,i]>(3*pi/2)){
     solution<-A1[j,]
     if(i!=1){
      a<-i-1
      b<-i
      a<-min(c(1:length(G1))[A1[j,]==A1[j,a]])
      b<-max(c(1:length(G1))[A1[j,]==A1[j,b]])
      solution[a:b]<-means[a,b]
      }
     if(i==1){
      b<-1;b<-max(c(1:length(G1))[A1[j,]==A1[j,b]])
      solution[1:b]<-0}
      A1<-rbind(A1,solution)
      }
     if(A1[j,i]>=pi&&A1[j,i]<=(3*pi/2)){      
      b<-i
      a<-i
      a<-max(c(1:length(G1))[A1[j,]==A1[j,a]])
      soli<-A1[j,i]
      while(soli>(pi/2)){
       b<-b-1
       if(b==0){soli<-0}
       if(b!=0){
        b<-min(c(1:length(G1))[A1[j,]==A1[j,b]])
        soli<-means[b,a]
    }      
       } # close while
      if(b==0){b<-1}
      solution<-A1[j,]
      solution[b:a]<-soli
      A1<-rbind(A1,solution) 
      b<-i
      a<-i 
      a<-min(c(1:length(G1))[A1[j,]==A1[j,a]])
      soli<-A1[j,i]
      exit<-0
      while(soli>(pi/2)&&exit==0){
       b<-b+1
       if(b>=(s+1)){
        exit<-1
        }
       if(b<=s){
        b<-max(c(1:length(G1))[A1[j,]==A1[j,b]])
        soli<-means[a,b]
        }
       } # close while
       if(exit==0){
        solution<-A1[j,]
        solution[a:b]<-soli
        A1<-rbind(A1,solution)
       } 
       if(i!=s){
        if(i==1){
         b<-2
         b<-max(c(1:length(G1))[A1[j,]==A1[j,b]])
         solution<-A1[j,]; solution[1:b]<-0
         A1<-rbind(A1,solution)
         } # close if i==1
        if(i!=1){
         a<-i-1
         a<-min(c(1:length(G1))[A1[j,]==A1[j,a]])
         b<-i+1
         b<-max(c(1:length(G1))[A1[j,]==A1[j,b]])
         solution<-A1[j,]
         solution[a:b]<-means[a,b]
         A1<-rbind(A1,solution)
         } # close if i!=1
        } # close if i!=s
      }# close if    
    }#close for(j) 
   if(nrow(A1)>fin){
    A1<-rbind(A1[(fin+1):nrow(A1),])
    }
   A1<-unique(A1,MARGIN=1) 
   }#close for(i)
  }#close if s!=0
 lA1[[s+1]]<-A1
 } # close for(s)

### G2
for(s in 0:q){
 lA2[[s+1]]<-list()
 for(p in (s+1):(q+1)){
  G2<-array()
  A2<-array()
  if(s<(p-1)&&p>1&&s<q&&p<(q+2)){
   G2<-data[(s+1):(p-1)]
   A2<-rbind(G2)
   for(i in 1:length(G2)){
    fin<-nrow(A2)
    for(j in 1:fin){
     if(A2[j,i]>(pi/2)&&A2[j,i]<(3*pi/2)){
      solution<-A2[j,]
      A2<-rbind(A2,solution)
      }# close if
     if(A2[j,i]>(3*pi/2)||A2[j,i]<(pi/2)){
      b<-i
      a<-i
      a<-max(c(1:length(G2))[A2[j,]==A2[j,a]])
      soli<-A2[j,i]
      exit<-0
      while((soli>(3*pi/2)||soli<(pi/2))&&exit==0){
       b<-b-1
       if(b==0){exit<-1}
       if(b>=1){
        b<-min(c(1:length(G2))[A2[j,]==A2[j,b]])
        soli<-means[(b+s),(a+s)]
        }
       }# close while
      if(exit==0){
       solution<-A2[j,]
       solution[b:a]<-soli
       A2<-rbind(A2,solution)
       } 
      b<-i
      a<-i 
      a<-min(c(1:length(G2))[A2[j,]==A2[j,a]])
      soli<-A2[j,i]
      exit<-0
      while((soli>(3*pi/2)||soli<(pi/2))&&exit==0){
       b<-b+1
       if(b>length(G2)){exit<-1}
       if(b<=length(G2)){
        b<-max(c(1:length(G2))[A2[j,]==A2[j,b]])
        soli<-means[(a+s),(b+s)]}      
       }# close while
      if(b==(length(G2)+1)){b<-length(G2)}
      if(exit==0){
       solution<-A2[j,]
       solution[a:b]<-soli
       A2<-rbind(A2,solution)
       } 
      if(i!=1&&i!=length(G2)&&i!=q){
       a<-i-1
       a<-min(c(1:length(G2))[A2[j,]==A2[j,a]])
       b<-i+1
       b<-max(c(1:length(G2))[A2[j,]==A2[j,b]])
       solution<-A2[j,]
       solution[a:b]<-means[(a+s),(b+s)]
       A2<-rbind(A2,solution)
       }# close if
      }# close if
     }#close for(j) 
     if(nrow(A2)>fin){A2<-rbind(A2[(fin+1):nrow(A2),])}
     A2<-unique(A2,MARGIN=1) 
    }#close for(i)
   }#close if G2
   lA2[[s+1]][[p]]<-A2
  }#close for(p)
 }#close for(s) 

### G3
for(p in 1:(q+1)){
 G3<-array()
 A3<-array()
 if(p!=(q+1)){
  G3<-data[p:q]
  A3<-rbind(G3)
  for(i in 1:length(G3)){
   fin<-nrow(A3)
   for(j in 1:fin){
    if(A3[j,i]>pi){
     solution<-A3[j,]
     A3<-rbind(A3,solution)
     }# close if
    if(A3[j,i]<(pi/2)){
     solution<-A3[j,]
     if(i!=q&& i!=length(G3)){
      a<-i
      b<-i+1
      a<-min(c(1:length(G3))[A3[j,]==A3[j,a]])
      b<-max(c(1:length(G3))[A3[j,]==A3[j,b]])
      solution[a:b]<-means[(a+(p-1)),(b+(p-1))]
      }
     if(i==length(G3)){
      b<-i;b<-min(c(1:length(G3))[A3[j,]==A3[j,b]])
      solution[b:length(G3)]<-2*pi
      }
     A3<-rbind(A3,solution)
     }# close if
    if(A3[j,i]>(pi/2)&&A3[j,i]<pi){
     b<-i
     a<-i
     a<-max(c(1:length(G3))[A3[j,]==A3[j,a]])
     soli<-A3[j,i]
     exit<-0
     while(soli<(3*pi/2)&&exit==0){
      b<-b-1
      if(b==0){exit<-1}
      if(b>=1){
       b<-min(c(1:length(G3))[A3[j,]==A3[j,b]])
       soli<-means[(b+(p-1)),(a+(p-1))]
       }
      }# close while
     if(exit==0){
      solution<-A3[j,]
      solution[b:a]<-soli
      A3<-rbind(A3,solution)
      } 
     b<-i
     a<-i
     a<-min(c(1:length(G3))[A3[j,]==A3[j,a]])
     soli<-A3[j,i]
     exit<-0
     while(soli<(3*pi/2)&&exit==0){
      b<-b+1
      if(b>length(G3)){exit<-1}
      if(b==(length(G3)+1)){
       soli<-2*pi
       exit<-0
       }
      if(b<=length(G3)&&exit==0){
       b<-max(c(1:length(G3))[A3[j,]==A3[j,b]])
       soli<-means[(a+(p-1)),(b+(p-1))]
       }
      }# close while
     if(b==(length(G3)+1)){b<-length(G3)}
     if(exit==0){
      solution<-A3[j,]
      solution[a:b]<-soli
      A3<-rbind(A3,solution)
      }
     if(i!=1&&i!=length(G3)){
      if(i==q){
       b<-q-1
       b<-min(c(1:length(G3))[A3[j,]==A3[j,b]])
       solution<-A3[j,]
       solution[b:q]<-2*pi
       A3<-rbind(A3,solution)
       }
      if(i!=q){
       a<-i-1
       a<-min(c(1:length(G3))[A3[j,]==A3[j,a]])
       b<-i+1
       b<-max(c(1:length(G3))[A3[j,]==A3[j,b]])
       solution<-A3[j,]; solution[a:b]<-means[(a+(p-1)),(b+(p-1))]
       A3<-rbind(A3,solution)
       }
      }# close if
     }# close if 
    }# close for(j) 
   if(nrow(A3)>fin){A3<-rbind(A3[(fin+1):nrow(A3),])}
   A3<-unique(A3,MARGIN=1) 
   }# close for(i)
  }# close if G3
 lA3[[p]]<-A3
 } # close for(p)
 
# applying PAVA 
lA1.PAVA<-list()
SCE1<-array()
for(i in 1:length(lA1)){
 if(!is.na(lA1[[i]][1])){ 
  lA1.PAVA[[i]]<-cirPAVA(lA1[[i]][1,])$data 
  SCE1[i]<-cirSCE(data[1:(i-1)],lA1.PAVA[[i]],mrl[1:(i-1)])
  if(nrow(lA1[[i]])>1){     
   for(n in 1:nrow(lA1[[i]])){    
    result1<-cirPAVA(lA1[[i]][n,])$data
    auxSCE<-cirSCE(data[1:(i-1)],result1,mrl[1:(i-1)])
    if(result1[length(result1)]<(pi/2)){
     if(auxSCE<SCE1[i]){
      lA1.PAVA[[i]]<-result1
      SCE1[i]<-auxSCE
      }# close if<SCE
     }# close if condition (5)
    }# close for solutions
   }# close if solutions
  if(lA1.PAVA[[i]][length(lA1.PAVA[[i]])]>(pi/2)){lA1.PAVA[[i]]<-rep(10000,length(lA1.PAVA[[i]]))}
  }# close check NA
 if(is.na(lA1[[i]][1])){lA1.PAVA[[i]]<-NA}
 }# close for(lA1) 
lA2.PAVA<-list()
SCE2<-matrix(ncol=length(lA2), nrow=length(lA2))
for(i in 1:length(lA2)){
 lA2.PAVA[[i]]<-list()
 for(j in i:length(lA2[[i]])){
  if(!is.na(lA2[[i]][[j]][1])){
   means2<-means[i:(j-1),i:(j-1)]
   lA2.PAVA[[i]][[j]]<-cirPAVA(lA2[[i]][[j]][1,])$data
   SCE2[i,j]<-cirSCE(data[i:(j-1)],lA2.PAVA[[i]][[j]],mrl[i:(j-1)])
   if(nrow(lA2[[i]][[j]])>1){
    for(n in 1:nrow(lA2[[i]][[j]])){
     result2<-cirPAVA(lA2[[i]][[j]][n,])$data
     auxSCE<-cirSCE(data[i:(j-1)],result2,mrl[i:(j-1)])
     if(result2[1]>(pi/2)&&result2[length(result2)]<(3*pi/2)){
      if(auxSCE<SCE2[i,j]){
       lA2.PAVA[[i]][[j]]<-result2
       SCE2[i,j]<-auxSCE
       }# close if<SCE
      }# close if conticion (5)
     }# close for solutions
    }# close if solutions
   if(lA2.PAVA[[i]][[j]][length(lA2.PAVA[[i]][[j]])]>(3*pi/2)||lA2.PAVA[[i]][[j]][1]<(pi/2)){lA2.PAVA[[i]][[j]]<-rep(10000,length(lA2.PAVA[[i]][[j]]))}
   }# close check NA
  if(is.na(lA2[[i]][[j]][1])){lA2.PAVA[[i]][[j]]<-NA}
  }# close for(j)
 }# close for(i)
lA3.PAVA<-list()
SCE3<-array()
for(i in 1:length(lA3)){
 if(!is.na(lA3[[i]][1])){
  means3<-means[i:nrow(means),i:ncol(means)]
  lA3.PAVA[[i]]<-cirPAVA(lA3[[i]][1,])$data
  SCE3[i]<-cirSCE(data[i:q],lA3.PAVA[[i]],mrl[i:q])
  if(nrow(lA3[[i]])>1){
   for(n in 1:nrow(lA3[[i]])){
    result3<-cirPAVA(lA3[[i]][n,])$data
    auxSCE<-cirSCE(data[i:q],result3,mrl[i:q])
    if(result3[1]>(3*pi/2)){
     if(auxSCE<SCE3[i]){
      lA3.PAVA[[i]]<-result3
      SCE3[i]<-auxSCE
      }# close if<SCE
     }# close if conticion (5)
    }# close for (n)
   }# close if solutions
  if(lA3.PAVA[[i]][1]<(3*pi/2)){lA3.PAVA[[i]]<-rep(10000,length(lA3.PAVA[[i]]))}
  }# close check NA
 if(is.na(lA3[[i]][1])){lA3.PAVA[[i]]<-NA}
 }# close for(i)

fis.tilde<-rbind(c(data,0))
for(i in 1:length(lA1)){
 for(j in i:length(lA1)){
  jointed<-c(lA1.PAVA[[i]],lA2.PAVA[[i]][[j]],lA3.PAVA[[j]])
  jointed<-jointed[!is.na(jointed)]
  SCE<-cirSCE(data,jointed,mrl)
  solution<-c(jointed,SCE)
  fis.tilde<-rbind(fis.tilde,solution)
  }# close for(j)
 }# close for(i) 

fis.tilde<-fis.tilde[2:nrow(fis.tilde),]
addition<-apply(fis.tilde,1,sum)
fis.tilde<-fis.tilde[addition<10000,]
finsol<-fis.tilde[fis.tilde[,ncol(fis.tilde)]==min(fis.tilde[,ncol(fis.tilde)]),]
CIRE<-finsol[1:(length(finsol)-1)]
SCE<-finsol[length(finsol)]
return(list(CIRE=CIRE,SCE=SCE))
}# end function CIRE

####################

means<-apply(data,1,cirmean)
resultant<-mrl(data)
ordmeans<-levels
names(ordmeans)<-means
sortmeans<-as.numeric(names(sort(ordmeans)))
ordresult<-levels
names(ordresult)<-resultant
sortresult<-as.numeric(names(sort(ordresult)))
n<-as.vector(table(levels))


k<-length(n)
q<-length(sortmeans)
listdata<-list(1)
listdata[[1]]<-sortmeans[1:n[1]]
  if( k > 1){
   for (l in 2:k){
    listdata[[l]]<-sortmeans[(cumsum(n)[(l-1)]+1):(cumsum(n)[l])]
    }
   }


ldata<-listdata
mrl<-sortresult
num<-n
nr<-prod(factorial(num))
nc<-sum(num)
n<-factorial(num)



indexes<-matrix(ncol=nc,nrow=nr)
mdata<-matrix(ncol=nc,nrow=nr)
aux<-1
for(p in 1:length(ldata)){
 l<-length(ldata[[p]])
 if(l==1){
  mdata[,aux]<-rep(ldata[[p]],nr)
  indexes[,aux]<-rep(aux,nr)
  aux<-aux+1
  }# end if level of size 1
 if(l>1){
  iaux<-c(aux:(aux+l-1))
  permu<-permn(ldata[[p]])
  ipermu<-permn(iaux)
  lp<-length(permu)
  nrep<-nr/lp
  if(p==1){
   maux<-matrix(nrow=lp,ncol=l)
   miaux<-matrix(nrow=lp,ncol=l)
   for(d in 1:lp){
    maux[d,]<-permu[[d]]
    miaux[d,]<-ipermu[[d]]
    }
   for(c in 1:nrep){
    mdata[(n[1]*(c-1)+1):(c*n[1]),1:l]<-maux
    indexes[(n[1]*(c-1)+1):(c*n[1]),1:l]<-miaux
    }
   aux<-aux+l
   } # close if p=1
  if(p!=1){
   repe<-prod(n[1:(p-1)])
   maux<-matrix(nrow=repe,ncol=l)
   miaux<-matrix(nrow=repe,ncol=l)
   for(a in 1:(nr/repe)){
    aux1<-a%%n[p]
    if(aux1==0){aux1<-n[p]}
    repe<-prod(n[1:(p-1)])
    maux<-matrix(nrow=repe,ncol=l)
    miaux<-matrix(nrow=repe,ncol=l)
    for(d in 1:repe){
     maux[d,]<-permu[[aux1]]
     miaux[d,]<-ipermu[[aux1]]
     }
    if(a>n[p]){aux1<-aux1+1}
    mdata[(((a-1)*repe)+1):(a*repe),aux:(aux+l-1)]<-maux
    indexes[(((a-1)*repe)+1):(a*repe),aux:(aux+l-1)]<-miaux
    } # close for a
   aux<-aux+l
   } # close if p!=1
  } # close if l>1
 } # close p

CIRE<-matrix(ncol=ncol(mdata),nrow=nrow(mdata))
SCE<-array()

for(m in 1:nrow(mdata)){
 data<-mdata[m,]
 ind<-indexes[m,]
 means<-mmeans(data)
 CIRE[m,]<-funCIRE(data,means,mrl)$CIRE
 SCE[m]<-funCIRE(data,means,mrl)$SCE
 if(isotropic==TRUE){
  aux<-1
  for (i in 2:length(data)){
   data<-c(data[2:length(data)],data[1])
   ind<-c(ind[2:length(ind)],ind[1])
   means<-mmeans(data)
   if(funCIRE(data,means,mrl)$SCE<SCE[m]){
    CIRE[m,]<-funCIRE(data,means,mrl)$CIRE
    SCE[m]<-funCIRE(data,means,mrl)$SCE
    indexes[m,]<-ind
    aux<-i
    } # close if<SCE
   } # close for(i)
  } # close isotropic=TRUE
 }# close for(m)

CIRESfinals<-CIRE[c(1:length(SCE))[SCE==min(SCE)],]
if(is.matrix(CIRESfinals)){eCIRE<-as.vector(unique(CIRESfinals))}
if(!is.matrix(CIRESfinals)){eCIRE<-CIRESfinals}
minSCE<-unique(SCE[SCE==min(SCE)])
finalindexes<-indexes[c(1:length(SCE))[SCE==min(SCE)],]
if(is.matrix(finalindexes)){names(eCIRE)<-finalindexes[1,]}
if(!is.matrix(finalindexes)){names(eCIRE)<-finalindexes}

sorted<-as.numeric(names(eCIRE))
names(sorted)<-eCIRE
sorted<-as.numeric(names(sort(sorted)))
cumform<-cumsum(num)
lCIRE<-list()
length(lCIRE)<-length(cumform)
lCIRE[[1]]<-sorted[1:cumform[1]]
if(length(cumform)>1){
 for(i in 2:length(cumform)){
  lCIRE[[i]]<-sorted[((cumform[(i-1)])+1):cumform[i]] 
  }
 }

lexit<-list(cirmeans=ldata,SCE=minSCE,CIRE=lCIRE)

if(graphic==TRUE){
 x11()
 plotcircularm(as.circular(sorted),stack=TRUE,bins=300,main="Circular Isotonic Regression Estimator",col=2,pch=c(1:length(sorted)))
 x11()
 plotcircularm(as.circular(mdata[1,]),stack=TRUE,bins=300,main="Unrestricted Estimator",col=1,pch=c(1:length(mdata[1,])))
 } # close if plot
return(lexit)

} ###### end function CIREi


