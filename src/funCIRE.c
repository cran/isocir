#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


double max(SEXP data){
	int ii;
	double max=0;
	for(ii=0;ii<length(data);ii++){
		max= REAL(data)[ii] > max ? REAL(data)[ii] : max;
	}
	return max;
}

double min(SEXP data){
	int ii;
	double min=2*M_PI;
	for(ii=0;ii<length(data);ii++){
		min= REAL(data)[ii] < min ? REAL(data)[ii] : min;
	}
	return min;
}

double sce(SEXP data1, SEXP data2, SEXP mrl){
	double aux=0;
	int ii;
	for(ii=0;ii<length(data1);ii++)aux=aux+(REAL(mrl)[ii]*(1-cos(REAL(data1)[ii]-REAL(data2)[ii])));
	return(aux);
}

double cirmean(SEXP data){
	int i;
	double A=sin(REAL(data)[0]), B=cos(REAL(data)[0]), result=0;	
	if(length(data)>1){
		for(i=1;i<length(data);i++){
			A=A+sin(REAL(data)[i]);
			B=B+cos(REAL(data)[i]);
		}
	}
	if ((A >= 0) && (B > 0))result=atan(A/B);
    	if (B < 0)result=atan(A/B) + M_PI;
	if ((A < 0) && (B > 0))result=atan(A/B) + (2 * M_PI);
    	return(result);
}

			


SEXP cirPAVA(SEXP datas, SEXP A, SEXP means, SEXP mrl, int group){
        
	SEXP aux, problems, aux2, aux3, result;
	int nProtected=0, condition5=0, len=length(datas);

	PROTECT(aux = allocVector(REALSXP,  len));
	++nProtected;
	PROTECT(result = allocVector(REALSXP,  len));
	++nProtected;

	PROTECT(problems = allocVector(REALSXP,  len-1));
	++nProtected;
	PROTECT(aux2 = allocVector(REALSXP,  len-1));
	++nProtected;

	PROTECT_INDEX ipx1;
	PROTECT_WITH_INDEX(aux3 = allocVector(REALSXP,  0), &ipx1);
     	++nProtected;
        

	int nrow, ii, i, j, sumproblems, z, a, b;
	double check, Sumce=1000000, nSumce, condition;
	
	nrow=length(A)/len;
	for(j=0;j<len;j++)REAL(result)[j]=10000;

	for(i=0;i<nrow;i++){
		
		for(j=0;j<len;j++)REAL(aux)[j]=REAL(A)[i*len+j];
	
		/* aux es data */
        	sumproblems=0;
    		for(j=0;j<len-1;j++){
			condition=REAL(aux)[j]-REAL(aux)[j+1];
    			if(condition>0){
				REAL(problems)[j]=1;
				REAL(aux2)[sumproblems]=j;
				sumproblems++;
    			}
    			if(condition<=0)REAL(problems)[j]=0;
    		}
	
		REPROTECT(aux3 = allocVector(REALSXP,  sumproblems), ipx1);
		for(j=0;j<length(aux3);j++)REAL(aux3)[j]=REAL(aux2)[j];
			
	
		/* aux2 y aux 3 son pos de R */
		while(sumproblems!=0){
			for(ii=0;ii<length(aux3);ii++){
				a=REAL(aux3)[ii];
				b=a+1;
				z=0;	
				check=REAL(aux)[a];
				while(check!=REAL(aux)[z])z++;
				a=z;

				z=length(aux)-1;
				check=REAL(aux)[b];
				while(check!=REAL(aux)[z])z--;
				b=z;

			    	for(j=a;j<=b;j++)REAL(aux)[j]=REAL(means)[a*len+b];
			/*	Rprintf("means ab: %f, ",REAL(means)[a*len+b]);*/
		
			}
			sumproblems=0;
			for(j=0;j<len-1;j++){
				condition=REAL(aux)[j]-REAL(aux)[j+1];
				if(condition>0){
					REAL(aux2)[sumproblems]=j;
					REAL(problems)[j]=1;
					sumproblems++;
				}
				if(condition<=0)REAL(problems)[j]=0;
			}
			REPROTECT(aux3 = allocVector(REALSXP,  sumproblems), ipx1);
			for(j=0;j<length(aux3);j++)REAL(aux3)[j]=REAL(aux2)[j];
		} 
		nSumce=sce(aux, datas, mrl);
	/*	Rprintf("SCE: %f",nSumce);*/
		if(group==1)condition5=REAL(aux)[length(aux)-1]<(M_PI/2);
		if(group==2)condition5=(REAL(aux)[0]>M_PI/2)&&(REAL(aux)[length(aux)-1]<(3*M_PI/2));
		if(group==3)condition5=REAL(aux)[0]>(3*M_PI/2);
		if((nSumce < Sumce && condition5)){
			Sumce= nSumce;
			for(j=0;j<len;j++)REAL(result)[j]=REAL(aux)[j];
		}
/*	Rprintf("\n");*/
	} 
      	UNPROTECT(nProtected);
     	return result;
}




SEXP funCIRE(SEXP data, SEXP means, SEXP mrl){
   R_len_t q=length(data);
   int i, j, jj, y, exit, s, p=1, a, b, z, ii, nProtected = 0, fin, d, f, count, nr, add;
       /*ncol, nrow, */
   double soli, check, aux=0.0;
   
    
   SEXP result;
   PROTECT(result = allocVector(VECSXP, 3));
   ++nProtected;

   SEXP lA1, lA2, lA2aux, lA3;
   SEXP solution;

 /*  Rprintf("data: ");
    for(ii=0;ii<q;ii++)Rprintf(" %f, ",REAL(data)[ii]);
    Rprintf("\n");*/

/* G1 */
   
   PROTECT(lA1 = allocVector(VECSXP, q));
   ++nProtected;

   SEXP A1, G1;

   PROTECT_INDEX ipx9;
   PROTECT_WITH_INDEX(G1 = allocVector(REALSXP,  0), &ipx9);
   ++nProtected;
  
   PROTECT_INDEX ipx1;
   PROTECT_WITH_INDEX(A1 = allocVector(REALSXP,  0), &ipx1);
   ++nProtected;

   PROTECT_INDEX ipx2;
   PROTECT_WITH_INDEX(solution = allocVector(REALSXP,  0), &ipx2);
   ++nProtected;

   SEXP auxA1;
   PROTECT_INDEX ipx3;
   PROTECT_WITH_INDEX(auxA1 = allocVector(REALSXP,  0), &ipx3);
   ++nProtected;

      for(s=0;s<q;s++){
   
    REPROTECT(A1 = allocVector(REALSXP,  s+1), ipx1);
    REPROTECT(G1 = allocVector(REALSXP,  s+1), ipx9);
    REPROTECT(solution = allocVector(REALSXP,  s+1), ipx2);
   

    for(ii=0;ii<length(A1);ii++)REAL(A1)[ii]=REAL(data)[ii];
    for(ii=0;ii<length(G1);ii++)REAL(G1)[ii]=REAL(data)[ii];

    for(i=s;i>=0;i--){
	  
	    fin=(length(A1)/(s+1));
	    for(j=0;j<fin;j++){
		    aux=REAL(A1)[j*(s+1)+i];
		 
		    if(aux<M_PI){
			    z=length(A1);
			    REPROTECT(auxA1 = allocVector(REALSXP,  length(A1)), ipx3);
			    for(ii=0;ii<length(A1);ii++)REAL(auxA1)[ii]=REAL(A1)[ii];
			    
			    REPROTECT(A1 = allocVector(REALSXP,  z+(s+1)), ipx1);
			    for(ii=0;ii<length(auxA1);ii++)REAL(A1)[ii]=REAL(auxA1)[ii];
			    for(ii=0;ii<s+1;ii++){
				    REAL(A1)[z]=REAL(auxA1)[j*(s+1)+ii];
				    z=z+1;
			    }
		    }		
		    if(aux<(2*M_PI) && aux>(3*M_PI/2)){
			    REPROTECT(solution = allocVector(REALSXP,  s+1), ipx2);
			    for(ii=0;ii<(s+1);ii++)REAL(solution)[ii]=REAL(A1)[j*(s+1)+ii];
			    if(i!=0){
				    a=i-1;
				    b=i;
				    
				    z=0;	     				    
				    check=REAL(A1)[(j*(s+1))+a];
				    while(check!=REAL(A1)[(j*(s+1))+z])z++;
				    a=z;					    

				   
				    z=(s+1)-1;
				    check=REAL(A1)[(j*(s+1))+b];
				    while(check!=REAL(A1)[(j*(s+1))+z])z--;
				    b=z;

				    for(ii=a;ii<=b;ii++)REAL(solution)[ii]=REAL(means)[a*length(data)+b];
			    }
			    if(i==0){
				    b=0;
				    
				    z=(s+1)-1;
				    check=REAL(A1)[(j*(s+1))+b];
				    while(check!=REAL(A1)[(j*(s+1))+z])z--;
				    b=z;
				    
				    for(ii=0;ii<=b;ii++)REAL(solution)[ii]=0.0;
			    }
			    add=1;
			    nr=(length(A1)-(fin*(s+1)))/(s+1);
			    if(length(A1)>fin*(s+1)){
				    for(d=0;d<nr;d++){
					    count=0;
					    for(f=0;f<(s+1);f++){
						    if(REAL(A1)[((d+fin)*(s+1))+f]==REAL(solution)[f])count++;
					    }
					    if(count==(s+1)){
						   add=0;
						   break;
					    }
				    }
			    }
		   	    if(add==1){
				    z=length(A1);
				    REPROTECT(auxA1 = allocVector(REALSXP,  length(A1)), ipx3);
				    for(ii=0;ii<length(A1);ii++)REAL(auxA1)[ii]=REAL(A1)[ii];
				    
				    REPROTECT(A1 = allocVector(REALSXP,  z+(s+1)), ipx1);
				    for(ii=0;ii<length(auxA1);ii++)REAL(A1)[ii]=REAL(auxA1)[ii];
				    for(ii=0;ii<length(solution);ii++){
					    REAL(A1)[z]=REAL(solution)[ii];
					    z=z+1;
				    }
			    }	    			    
		    }
		    if(aux>=M_PI && aux<=(3*M_PI/2)){
			    b=i;
			    a=i;
			    
			    z=(s+1)-1;
			    check=REAL(A1)[(j*(s+1))+a];
			    while(check!=REAL(A1)[(j*(s+1))+z])z--;
			    a=z;

			    soli=REAL(A1)[(j*(s+1))+i];

 			    while(soli>(M_PI/2)){
				    b=b-1;
				    if(b<=-1) soli=0;
				    if(b>=0){					    
					    z=0;
					    check=REAL(A1)[(j*(s+1))+b];
					    while(check!=REAL(A1)[(j*(s+1))+z])z++;
					    b=z;

					    soli=REAL(means)[b*length(data)+a];
				    }
			    }

			    
			    if(b<=-1) b=0;
			    REPROTECT(solution = allocVector(REALSXP,  s+1), ipx2);
			    for(ii=0;ii<(s+1);ii++)REAL(solution)[ii]=REAL(A1)[j*(s+1)+ii];
			   
			    for(ii=b;ii<=a;ii++)REAL(solution)[ii]=soli;

			    add=1;
			    nr=(length(A1)-(fin*(s+1)))/(s+1);
			    if(length(A1)>fin*(s+1)){
				    for(d=0;d<nr;d++){
					    count=0;
					    for(f=0;f<(s+1);f++){
						    if(REAL(A1)[((d+fin)*(s+1))+f]==REAL(solution)[f])count++;
					    }
					    if(count==(s+1)){
						   add=0;
						   break;
					    }
				    }
			    }
		   	    if(add==1){
				    z=length(A1);			    
				    REPROTECT(auxA1 = allocVector(REALSXP,  length(A1)), ipx3);	   
				    for(ii=0;ii<length(A1);ii++)REAL(auxA1)[ii]=REAL(A1)[ii];
				    REPROTECT(A1 = allocVector(REALSXP,  z+(s+1)), ipx1);
				    for(ii=0;ii<length(auxA1);ii++)REAL(A1)[ii]=REAL(auxA1)[ii];
				    for(ii=0;ii<length(solution);ii++){
					    REAL(A1)[z]=REAL(solution)[ii];
					    z=z+1;
				    }
			    }
			    b=i;
			    a=i;
			   
			    z=0;
			    check=REAL(A1)[(j*(s+1))+a];
			    while(check!=REAL(A1)[(j*(s+1))+z])z++;
			    a=z;

			    soli=REAL(A1)[(j*(s+1))+i];
			    exit=0;
			    while(soli>M_PI/2 && exit==0){
 				    b=b+1;
  				    if(b>=s+1)exit=1;
   				    if(b<=s){   					   
   					    z=(s+1)-1;
					    check=REAL(A1)[(j*(s+1))+b];
					    while(check!=REAL(A1)[(j*(s+1))+z])z--;
					    b=z;

   					    soli=REAL(means)[a*length(data)+b];
   				    }
			    }
			    if(exit==0){
				    for(ii=0;ii<(s+1);ii++)REAL(solution)[ii]=REAL(A1)[j*(s+1)+ii];
				    for(ii=a;ii<=b;ii++)REAL(solution)[ii]=soli;


				     add=1;
	 			     nr=(length(A1)-(fin*(s+1)))/(s+1);
	 			     if(length(A1)>fin*(s+1)){
	 				     for(d=0;d<nr;d++){
	 					     count=0;
	 					     for(f=0;f<(s+1);f++){
	 						     if(REAL(A1)[((d+fin)*(s+1))+f]==REAL(solution)[f])count++;
	 					     }
	 					     if(count==(s+1)){
	  						     add=0;
	  						     break;
	 					     }
	 				     }
	 			     }
	 			     if(add==1){
	 				     z=length(A1);
	 				     REPROTECT(auxA1 = allocVector(REALSXP,  length(A1)), ipx3);
	 				     for(ii=0;ii<length(A1);ii++)REAL(auxA1)[ii]=REAL(A1)[ii];
	 				     
	 				     REPROTECT(A1 = allocVector(REALSXP,  z+(s+1)), ipx1);
	 				     for(ii=0;ii<length(auxA1);ii++)REAL(A1)[ii]=REAL(auxA1)[ii];
	 				     for(ii=0;ii<length(solution);ii++){
	 					     REAL(A1)[z]=REAL(solution)[ii];
	 					     z=z+1;
	 				     }
				     }
			    }
			    if(i!=s){
				    if(i==0){
					    b=1;
					    z=(s+1)-1;
					    check=REAL(A1)[(j*(s+1))+b];
					    while(check!=REAL(A1)[(j*(s+1))+z])z--;
					    b=z;

					    for(ii=0;ii<(s+1);ii++)REAL(solution)[ii]=REAL(A1)[j*(s+1)+ii];
					    for(ii=0;ii<=b;ii++)REAL(solution)[ii]=0;

					     add=1;
		 			     nr=(length(A1)-(fin*(s+1)))/(s+1);
		 			     if(length(A1)>fin*(s+1)){
		 				     for(d=0;d<nr;d++){
		 					     count=0;
		 					     for(f=0;f<(s+1);f++){
		 						     if(REAL(A1)[((d+fin)*(s+1))+f]==REAL(solution)[f])count++;
		 					     }
		 					     if(count==(s+1)){
		  						     add=0;
		  						     break;
		 					     }
		 				     }
		 			     }
		 			     if(add==1){
	 					     z=length(A1);
	 					     REPROTECT(auxA1 = allocVector(REALSXP,  length(A1)), ipx3);
	 					     for(ii=0;ii<length(A1);ii++)REAL(auxA1)[ii]=REAL(A1)[ii];
	 					     REPROTECT(A1 = allocVector(REALSXP,  z+(s+1)), ipx1);
	 					     for(ii=0;ii<length(auxA1);ii++)REAL(A1)[ii]=REAL(auxA1)[ii];
	 					     for(ii=0;ii<length(solution);ii++){
	 						     REAL(A1)[z]=REAL(solution)[ii];
	 						     z=z+1;
	 					     }
					     }
				    }
				    if(i!=0){
					    a=i-1;
					    z=0;
 					    check=REAL(A1)[(j*(s+1))+a];
					    while(check!=REAL(A1)[(j*(s+1))+z])z++;
					    a=z;

					    b=i+1;					   
					    z=(s+1)-1;
					    check=REAL(A1)[(j*(s+1))+b];
					    while(check!=REAL(A1)[(j*(s+1))+z])z--;
					    b=z;

					    for(ii=0;ii<(s+1);ii++)REAL(solution)[ii]=REAL(A1)[j*(s+1)+ii];
       					    for(ii=a;ii<=b;ii++)REAL(solution)[ii]=REAL(means)[a*length(data)+b];

					     add=1;
		 			     nr=(length(A1)-(fin*(s+1)))/(s+1);
		 			     if(length(A1)>fin*(s+1)){
		 				     for(d=0;d<nr;d++){
		 					     count=0;
		 					     for(f=0;f<(s+1);f++){
		 						     if(REAL(A1)[((d+fin)*(s+1))+f]==REAL(solution)[f])count++;
		 					     }
		 					     if(count==(s+1)){
		  						     add=0;
		  						     break;
		 					     }
		 				     }
		 			     }
		 			     if(add==1){
	 					     z=length(A1);
	 					     REPROTECT(auxA1 = allocVector(REALSXP,  length(A1)), ipx3);
	 					     for(ii=0;ii<length(A1);ii++)REAL(auxA1)[ii]=REAL(A1)[ii];
	 					     REPROTECT(A1 = allocVector(REALSXP,  z+(s+1)), ipx1);
	 					     for(ii=0;ii<length(auxA1);ii++)REAL(A1)[ii]=REAL(auxA1)[ii];
	 					     for(ii=0;ii<length(solution);ii++){
	 						     REAL(A1)[z]=REAL(solution)[ii];
	 						     z=z+1;
	 					     } 
					     }
				    }
			    }
			   
		    }
	    }
	    if((length(A1)/(s+1))>fin){
		     REPROTECT(auxA1 = allocVector(REALSXP,  length(A1)-(fin*(s+1))), ipx3);
		     z=0;
		     for(ii=fin*(s+1);ii<length(A1);ii++){
			     REAL(auxA1)[z]=REAL(A1)[ii];
			     z=z+1;
		     }
		     REPROTECT(A1 = allocVector(REALSXP,  length(auxA1)), ipx1);
		     for(ii=0;ii<length(auxA1);ii++)REAL(A1)[ii]=REAL(auxA1)[ii];
	    }
	  
    } 
    
   /* SEXP m1;
    int nrow, ncol;
    ncol=s+1;
    nrow=length(A1)/(s+1);
    PROTECT(m1 = allocMatrix(REALSXP, nrow, ncol));
    ++nProtected;
    for(ii=0;ii<nrow;ii++){
	    for(jj=0;jj<ncol;jj++){
		    REAL(m1)[ii+(nrow*jj)]=REAL(A1)[ii*ncol+jj];
	    }
    }*/ 
    /* fin primera parte */

    /* cirPAVA */

    SEXP means1;
    PROTECT(means1 = allocMatrix(REALSXP, s+1, s+1));
    ++nProtected;
    for(ii=0;ii<s+1;ii++){
	    for(jj=0;jj<s+1;jj++){
		    REAL(means1)[ii+(s+1)*jj]=REAL(means)[ii+length(data)*jj];
	    }
    } 
   solution=cirPAVA(G1, A1, means1, mrl, 1);
               
    SET_VECTOR_ELT(lA1, s, solution);

   } /* fin for s */



/* G2 */

   PROTECT(lA2 = allocVector(VECSXP, q));
   ++nProtected;
   SEXP A2, auxA2, G2;

   PROTECT_INDEX ipx8;
   PROTECT_WITH_INDEX(lA2aux = allocVector(VECSXP, q), &ipx8);
   ++nProtected;

   PROTECT_INDEX ipx4;
   PROTECT_WITH_INDEX(A2 = allocVector(REALSXP,  0), &ipx4);
   ++nProtected;

   PROTECT_INDEX ipx5;
   PROTECT_WITH_INDEX(auxA2 = allocVector(REALSXP,  0), &ipx5);
   ++nProtected;

    PROTECT_INDEX ipx10;
   PROTECT_WITH_INDEX(G2 = allocVector(REALSXP,  0), &ipx10);
   ++nProtected;

  
   for(s=0;s<q;s++){   
   
     REPROTECT(lA2aux = allocVector(VECSXP,  q), ipx8);
    for(p=s+1;p<=q;p++){
	   
	   
	     REPROTECT(A2 = allocVector(REALSXP,  p-s), ipx4);
	     REPROTECT(G2 = allocVector(REALSXP,  p-s), ipx10);
	     REPROTECT(solution = allocVector(REALSXP,  p-s), ipx2);

	    

	     z=0;
	     for(ii=s;ii<p;ii++){
		     REAL(A2)[z]=REAL(data)[ii];
		     z=z+1;
	     }

	      z=0;
	     for(ii=s;ii<p;ii++){
		     REAL(G2)[z]=REAL(data)[ii];
		     z=z+1;
	     }
	     for(i=0;i<(p-s);i++){		  
		     fin=length(A2)/(p-s);
		     for(j=0;j<fin;j++){			    
			     aux=REAL(A2)[j*(p-s)+i];			    		     
			     if(aux>M_PI/2 && aux<3*M_PI/2){
				     for(ii=0;ii<(p-s);ii++)REAL(solution)[ii]=REAL(A2)[j*(p-s)+ii];

				     z=length(A2);
				     REPROTECT(auxA2 = allocVector(REALSXP,  length(A2)), ipx5);
				     for(ii=0;ii<length(A2);ii++)REAL(auxA2)[ii]=REAL(A2)[ii];
				     REPROTECT(A2 = allocVector(REALSXP,  z+(p-s)), ipx4);
				     for(ii=0;ii<length(auxA2);ii++)REAL(A2)[ii]=REAL(auxA2)[ii];
				     for(ii=0;ii<length(solution);ii++){
					     REAL(A2)[z]=REAL(solution)[ii];
					     z=z+1;
				     }

			     }
			     if(aux>3*M_PI/2 || aux<M_PI/2){
				     b=i;
				     a=i;
				     z=(p-s)-1;
				     check=REAL(A2)[(j*(p-s))+a];
				     while(check!=REAL(A2)[(j*(p-s))+z])z--;
				     a=z;
				     
				     soli=aux;
				     exit=0;
				    
				     while((soli>3*M_PI/2 || soli<M_PI/2) && exit==0){
					     b=b-1;
					     if(b<0)exit=1;
					     if(b>=0){
						     z=0;
						     check=REAL(A2)[(j*(p-s))+b];
						     while(check!=REAL(A2)[(j*(p-s))+z])z++;
						     b=z;
						     
						     soli=REAL(means)[(b+s)+length(data)*(a+s)];
					     }
				     }
				     
				     if(exit==0){
					     for(ii=0;ii<(p-s);ii++)REAL(solution)[ii]=REAL(A2)[j*(p-s)+ii];
					     for(ii=b;ii<=a;ii++)REAL(solution)[ii]=soli;

					      add=1;
		  			      nr=(length(A2)-(fin*(p-s)))/(p-s);
		  			      if(length(A2)>fin*(p-s)){
		  				      for(d=0;d<nr;d++){
		  					      count=0;
		  					      for(f=0;f<(p-s);f++){
		  						      if(REAL(A2)[((d+fin)*(p-s))+f]==REAL(solution)[f])count++;
		  					      }
		  					      if(count==(p-s)){
		   						      add=0;
		   						      break;
		  					      }
		  				      }
		  			      }
		  			      if(add==1){
	 					      z=length(A2);
	 					      REPROTECT(auxA2 = allocVector(REALSXP,  length(A2)), ipx5);
	 					      for(ii=0;ii<length(A2);ii++)REAL(auxA2)[ii]=REAL(A2)[ii];
	 					      REPROTECT(A2 = allocVector(REALSXP,  z+(p-s)), ipx4);
	 					      for(ii=0;ii<length(auxA2);ii++)REAL(A2)[ii]=REAL(auxA2)[ii];
	 					      for(ii=0;ii<length(solution);ii++){
	 						      REAL(A2)[z]=REAL(solution)[ii];
	 						      z=z+1;
	 					      }
					      }
				     } 
				     b=i;
				     a=i;
				     z=0;
				     check=REAL(A2)[(j*(p-s))+a];
				     while(check!=REAL(A2)[(j*(p-s))+z])z++;
				     a=z;
				     soli=aux;
				     exit=0;
				     while((soli>3*M_PI/2 || soli<M_PI/2) && exit==0){
					     b=b+1;
					     if(b>=(p-s))exit=1;
					     if(b<(p-s)){
						     z=(p-s)-1;
						     check=REAL(A2)[(j*(p-s))+b];
						     while(check!=REAL(A2)[(j*(p-s))+z])z--;
						     b=z;
						     
						     soli=REAL(means)[(a+s)+length(data)*(b+s)];
					     }
				     }
				     
				     if(b==(p-s))b=(p-s)-1;
				     if(exit==0){
					     for(ii=0;ii<(p-s);ii++)REAL(solution)[ii]=REAL(A2)[j*(p-s)+ii];
					     for(ii=a;ii<=b;ii++)REAL(solution)[ii]=soli;

					      add=1;
					      nr=(length(A2)-(fin*(p-s)))/(p-s);
		  			      if(length(A2)>fin*(p-s)){
		  				      for(d=0;d<nr;d++){
		  					      count=0;
		  					      for(f=0;f<(p-s);f++){
		  						      if(REAL(A2)[((d+fin)*(p-s))+f]==REAL(solution)[f])count++;
		  					      }
		  					      if(count==(p-s)){
		   						      add=0;
		   						      break;
		  					      }
		  				      }
		  			      }
		  			      if(add==1){					     
					     z=length(A2);
					     REPROTECT(auxA2 = allocVector(REALSXP,  length(A2)), ipx5);
					     for(ii=0;ii<length(A2);ii++)REAL(auxA2)[ii]=REAL(A2)[ii];
					     
					     REPROTECT(A2 = allocVector(REALSXP,  z+(p-s)), ipx4);
					     for(ii=0;ii<length(auxA2);ii++)REAL(A2)[ii]=REAL(auxA2)[ii];
					     for(ii=0;ii<length(solution);ii++){
						     REAL(A2)[z]=REAL(solution)[ii];
						     z=z+1;
					     } 
					      }
				     }
				     if((i!=0) && (i!=(p-s)-1) && (i!=q-1)){
					     a=i-1;
					     z=0;					     
					     check=REAL(A2)[(j*(p-s))+a];
					     while(check!=REAL(A2)[(j*(p-s))+z])z++;
					     a=z;
					     
					     b=i+1;
					     z=(p-s)-1;
					     check=REAL(A2)[(j*(p-s))+b];
					     while(check!=REAL(A2)[(j*(p-s))+z])z--;
					     b=z;
					    
					     for(ii=0;ii<(p-s);ii++)REAL(solution)[ii]=REAL(A2)[j*(p-s)+ii];
					     for(ii=a;ii<=b;ii++)REAL(solution)[ii]=REAL(means)[(a+s)+length(data)*(b+s)];	     
					    
					      add=1;
		  			      nr=(length(A2)-(fin*(p-s)))/(p-s);
		  			      if(length(A2)>fin*(p-s)){
		  				      for(d=0;d<nr;d++){
		  					      count=0;
		  					      for(f=0;f<(p-s);f++){
		  						      if(REAL(A2)[((d+fin)*(p-s))+f]==REAL(solution)[f])count++;
		  					      }
		  					      if(count==(p-s)){
		   						      add=0;
		   						      break;
		  					      }
		  				      }
		  			      }
		  			      if(add==1){				     				     
					     z=length(A2);
					     REPROTECT(auxA2 = allocVector(REALSXP,  length(A2)), ipx5);
					     for(ii=0;ii<length(A2);ii++)REAL(auxA2)[ii]=REAL(A2)[ii];
					     REPROTECT(A2 = allocVector(REALSXP,  z+(p-s)), ipx4);
					     for(ii=0;ii<length(auxA2);ii++)REAL(A2)[ii]=REAL(auxA2)[ii];
					     for(ii=0;ii<length(solution);ii++){
						     REAL(A2)[z]=REAL(solution)[ii];
						     z=z+1;
					     } 
					      }
				     }

			     } 
		     }
		    
		     if((length(A2)/(p-s))>fin){
			     REPROTECT(auxA2 = allocVector(REALSXP,  length(A2)-(fin*(p-s))), ipx5);
			     z=0;
			     for(ii=fin*(p-s);ii<length(A2);ii++){
				     REAL(auxA2)[z]=REAL(A2)[ii];
				     z=z+1;
			     }
			     REPROTECT(A2 = allocVector(REALSXP,  length(auxA2)), ipx4);
			     for(ii=0;ii<length(auxA2);ii++)REAL(A2)[ii]=REAL(auxA2)[ii];
		     } 
		
	     } 
     	   
	     /*SEXP m2; 
	     ncol=p-s;
	     nrow=length(A2)/(p-s);
	     PROTECT(m2 = allocMatrix(REALSXP, nrow, ncol));
	     ++nProtected;
	     for(ii=0;ii<nrow;ii++){
	 	     for(jj=0;jj<ncol;jj++){
			    
	 		     REAL(m2)[ii+(nrow*jj)]=REAL(A2)[ii*ncol+jj];
			   
	 	     }
	     } */
	     
	     SEXP mrl2;
	     PROTECT(mrl2 = allocVector(REALSXP, p-s));
	     ++nProtected;
	     z=0;
	     for(ii=s; ii<p;ii++){
		     REAL(mrl2)[z]=REAL(mrl)[ii];
		     z++;
	     } 

	     SEXP means2;
	     PROTECT(means2 = allocMatrix(REALSXP, p-s, p-s));
	     ++nProtected;
	    
	     y=0;
	     for(ii=s;ii<p;ii++){
		     z=0;
       		     for(jj=s;jj<p;jj++){
	 		     REAL(means2)[y+(p-s)*z]=REAL(means)[ii+length(data)*jj];
			     z++;
	 	     }
		     y++;
	     }  
	     solution=cirPAVA(G2, A2, means2, mrl2, 2);
	     

   	     SET_VECTOR_ELT(lA2aux, p-1, solution);	       
    } /*close p*/
    SET_VECTOR_ELT(lA2, s, lA2aux);
   } /* close s*/

   
 /* G3 */

   PROTECT(lA3 = allocVector(VECSXP, q));
   ++nProtected;
   SEXP A3, auxA3, G3;

   PROTECT_INDEX ipx6;
   PROTECT_WITH_INDEX(A3 = allocVector(REALSXP,  0), &ipx6);
   ++nProtected;

   PROTECT_INDEX ipx7;
   PROTECT_WITH_INDEX(auxA3 = allocVector(REALSXP,  0), &ipx7);
   ++nProtected;
   
   PROTECT_INDEX ipx11;
   PROTECT_WITH_INDEX(G3 = allocVector(REALSXP,  0), &ipx11);
   ++nProtected;

   for(p=0;p<q;p++){

    
    REPROTECT(A3 = allocVector(REALSXP, q-p), ipx6);
    i=0;
    for(ii=p;ii<q;ii++){
     REAL(A3)[i]=REAL(data)[ii];
     i=i+1;
    }
    REPROTECT(G3 = allocVector(REALSXP, q-p), ipx11);
    i=0;
    for(ii=p;ii<q;ii++){
     REAL(G3)[i]=REAL(data)[ii];
     i=i+1;
    }

    REPROTECT(solution = allocVector(REALSXP,  q-p), ipx2);
    for(i=0;i<(q-p);i++){
	    fin=length(A3)/(q-p);
	    for(j=0;j<fin;j++){
		    aux=REAL(A3)[j*(q-p)+i];
		    if(aux>M_PI){
			    z=length(A3);
		
			    REPROTECT(auxA3 = allocVector(REALSXP,  length(A3)), ipx7);
			    for(ii=0;ii<length(A3);ii++)REAL(auxA3)[ii]=REAL(A3)[ii];

                            REPROTECT(A3 = allocVector(REALSXP,  z+(q-p)), ipx6);
			    for(ii=0;ii<length(auxA3);ii++)REAL(A3)[ii]=REAL(auxA3)[ii];
			    for(ii=0;ii<(q-p);ii++){
				    REAL(A3)[z]=REAL(auxA3)[j*(q-p)+ii];
				    z=z+1;
			    }
		    }
		    if(aux<M_PI/2){
       			    REPROTECT(solution = allocVector(REALSXP,  q-p), ipx2);
			    for(ii=0;ii<(q-p);ii++)REAL(solution)[ii]=REAL(A3)[j*(q-p)+ii];

			    if(i!=(q-1) && i!=((q-p)-1)){
				    a=i;
				    b=i+1;

				    z=0;	     				    
				    check=REAL(A3)[(j*(q-p))+a];
				    while(check!=REAL(A3)[(j*(q-p))+z])z++;
				    a=z;

				    z=(q-p)-1;
				    check=REAL(A3)[(j*(q-p))+b];
				    while(check!=REAL(A3)[(j*(q-p))+z])z--;
				    b=z;

				    for(ii=a;ii<=b;ii++)REAL(solution)[ii]=REAL(means)[(a+p)+length(data)*(b+p)];
			    }
			    if(i==((q-p)-1)){
				    b=i;
				    z=0;	     				    
				    check=REAL(A3)[(j*(q-p))+b];
				    while(check!=REAL(A3)[(j*(q-p))+z])z++;
				    b=z;

				    for(ii=b;ii<length(solution);ii++)REAL(solution)[ii]=2*M_PI;
			    }
			    add=1;
			    nr=(length(A3)-(fin*(q-p)))/(q-p);
			    if(length(A3)>fin*(q-p)){
				    for(d=0;d<nr;d++){
					    count=0;
					    for(f=0;f<(q-p);f++){
						    if(REAL(A3)[((d+fin)*(q-p))+f]==REAL(solution)[f])count++;
					    }
					    if(count==(q-p)){
						    add=0;
						    break;
					    }
				    }
			    }
			    if(add==1){

			    z=length(A3);
			 
			    REPROTECT(auxA3 = allocVector(REALSXP,  length(A3)), ipx7);
			    for(ii=0;ii<length(A3);ii++)REAL(auxA3)[ii]=REAL(A3)[ii];

                            REPROTECT(A3 = allocVector(REALSXP,  z+(q-p)), ipx6);

			    for(ii=0;ii<length(auxA3);ii++)REAL(A3)[ii]=REAL(auxA3)[ii];
			    for(ii=0;ii<q-p;ii++){
				    REAL(A3)[z]=REAL(solution)[ii];
				    z=z+1;
			    }
			    }
		    }
		    if(aux>M_PI/2 && aux<M_PI){

			    REPROTECT(solution = allocVector(REALSXP,  q-p), ipx2);

			    b=i;
			    a=i;

			    z=(q-p)-1;
			    check=REAL(A3)[(j*(q-p))+a];
			    while(check!=REAL(A3)[(j*(q-p))+z])z--;
			    a=z;

			    soli=aux;
			    exit=0;
			    while(soli<(3*M_PI/2) && exit==0){
				    b=b-1;
				    if(b<0)exit=1;
				    if(b>=0){
					    z=0;
					    check=REAL(A3)[(j*(q-p))+b];
					    while(check!=REAL(A3)[(j*(q-p))+z])z++;
					    b=z;
	
					    soli=REAL(means)[(b+p)+length(data)*(a+p)];
				    }
			    }
			    if(exit==0){
				    
				    for(ii=0;ii<(q-p);ii++)REAL(solution)[ii]=REAL(A3)[j*(q-p)+ii];
				    for(ii=b;ii<=a;ii++)REAL(solution)[ii]=soli;

				    add=1;
				    nr=(length(A3)-(fin*(q-p)))/(q-p);
				    if(length(A3)>fin*(q-p)){
					    for(d=0;d<nr;d++){
						    count=0;
						    for(f=0;f<(q-p);f++){
							    if(REAL(A3)[((d+fin)*(q-p))+f]==REAL(solution)[f])count++;
						    }
						    if(count==(q-p)){
							    add=0;
							    break;
						    }
					    }
				    }
				    if(add==1){

				    z=length(A3);
				  
				    REPROTECT(auxA3 = allocVector(REALSXP,  length(A3)), ipx7);
				    for(ii=0;ii<length(A3);ii++)REAL(auxA3)[ii]=REAL(A3)[ii];
				    
	   			    REPROTECT(A3 = allocVector(REALSXP,  z+(q-p)), ipx6);
				    for(ii=0;ii<length(auxA3);ii++)REAL(A3)[ii]=REAL(auxA3)[ii];
				    for(ii=0;ii<(q-p);ii++){
					    REAL(A3)[z]=REAL(solution)[ii];
					    z=z+1;
				    }
				    }
			    }
			    b=i;
			    a=i;

			    z=0;
			    check=REAL(A3)[(j*(q-p))+a];
			    while(check!=REAL(A3)[(j*(q-p))+z])z++;
			    a=z;
			    
			    soli=aux;
			    exit=0;
			    while(soli<(3*M_PI/2) && exit==0){
				    b=b+1;
				    if(b>(q-p))exit=1;
				    if(b==(q-p)){
					    soli=2*M_PI;
					    exit=0;					
				    }
				    if(b<(q-p) && exit==0){
					    z=(q-p)-1;
					    check=REAL(A3)[(j*(q-p))+b];
					    while(check!=REAL(A3)[(j*(q-p))+z])z--;
					    b=z;
					    
					    soli=REAL(means)[(a+p)+length(data)*(b+p)];
				    }
			    }
			  
			    if(b==(q-p))b=b-1;
			    if(exit==0){
				   
				    for(ii=0;ii<(q-p);ii++)REAL(solution)[ii]=REAL(A3)[j*(q-p)+ii];
				    for(ii=a;ii<=b;ii++)REAL(solution)[ii]=soli;

				    add=1;
				    nr=(length(A3)-(fin*(q-p)))/(q-p);
				    if(length(A3)>fin*(q-p)){
					    for(d=0;d<nr;d++){
						    count=0;
						    for(f=0;f<(q-p);f++){
							    if(REAL(A3)[((d+fin)*(q-p))+f]==REAL(solution)[f])count++;
						    }
						    if(count==(q-p)){
							    add=0;
							    break;
						    }
					    }
				    }
				    if(add==1){
				   
				    z=length(A3);
				    REPROTECT(auxA3 = allocVector(REALSXP,  z), ipx7);
				    for(ii=0;ii<length(A3);ii++)REAL(auxA3)[ii]=REAL(A3)[ii];
				    
	   			    REPROTECT(A3 = allocVector(REALSXP,  z+(q-p)), ipx6);
				    for(ii=0;ii<length(auxA3);ii++)REAL(A3)[ii]=REAL(auxA3)[ii];
				    for(ii=0;ii<(q-p);ii++){
					    REAL(A3)[z]=REAL(solution)[ii];
					    z=z+1;
				    }
				    }
			    }
			    if(i!=0 && i!=((q-p)-1)){
				    if(i==q-1){
					    b=q-2;
					    z=0;
					    check=REAL(A3)[(j*(q-p))+b];
					    while(check!=REAL(A3)[(j*(q-p))+z])z++;
					    b=z;
					    
					    
					    for(ii=0;ii<(q-p);ii++)REAL(solution)[ii]=REAL(A3)[j*(q-p)+ii];
					    for(ii=b;ii<q;ii++)REAL(solution)[ii]=2*M_PI;

					    add=1;
					    nr=(length(A3)-(fin*(q-p)))/(q-p);
					    if(length(A3)>fin*(q-p)){
						    for(d=0;d<nr;d++){
							    count=0;
							    for(f=0;f<(q-p);f++){
								    if(REAL(A3)[((d+fin)*(q-p))+f]==REAL(solution)[f])count++;
							    }
							    if(count==(q-p)){
								    add=0;
								    break;
							    }
						    }
					    }
					    if(add==1){

					    z=length(A3);				
					    REPROTECT(auxA3 = allocVector(REALSXP,  length(A3)), ipx7);
					    for(ii=0;ii<length(A3);ii++)REAL(auxA3)[ii]=REAL(A3)[ii];
					    
					    REPROTECT(A3 = allocVector(REALSXP,  z+(q-p)), ipx6);
					    for(ii=0;ii<length(auxA3);ii++)REAL(A3)[ii]=REAL(auxA3)[ii];
					    for(ii=0;ii<(q-p);ii++){
						    REAL(A3)[z]=REAL(solution)[ii];
						    z=z+1;
					    }
					    }
				    }
				    if(i!=q-1){
					    a=i-1;
					    z=0;
					    check=REAL(A3)[(j*(q-p))+a];
					    while(check!=REAL(A3)[(j*(q-p))+z])z++;
					    a=z;

					    b=i+1;
					    z=(q-p)-1;
					    check=REAL(A3)[(j*(q-p))+b];
					    while(check!=REAL(A3)[(j*(q-p))+z])z--;
					    b=z;

					    for(ii=0;ii<(q-p);ii++)REAL(solution)[ii]=REAL(A3)[j*(q-p)+ii];
					    for(ii=a;ii<=b;ii++)REAL(solution)[ii]=REAL(means)[(a+p)+length(data)*(b+p)];

					    add=1;
					    nr=(length(A3)-(fin*(q-p)))/(q-p);
					    if(length(A3)>fin*(q-p)){
						    for(d=0;d<nr;d++){
							    count=0;
							    for(f=0;f<(q-p);f++){
								    if(REAL(A3)[((d+fin)*(q-p))+f]==REAL(solution)[f])count++;
							    }
							    if(count==(q-p)){
								    add=0;
								    break;
							    }
						    }
					    }
					    if(add==1){

					    z=length(A3);					 
					    REPROTECT(auxA3 = allocVector(REALSXP,  length(A3)), ipx7);
					    for(ii=0;ii<length(A3);ii++)REAL(auxA3)[ii]=REAL(A3)[ii];
					    
					    REPROTECT(A3 = allocVector(REALSXP,  z+(q-p)), ipx6);
					    for(ii=0;ii<length(auxA3);ii++)REAL(A3)[ii]=REAL(auxA3)[ii];
					    for(ii=0;ii<(q-p);ii++){
						    REAL(A3)[z]=REAL(solution)[ii];
						    z=z+1;
					    }
					    }
				    }
			    }			 
		    }
	    } 
	    if((length(A3)/(q-p))>fin){
		     REPROTECT(auxA3 = allocVector(REALSXP,  length(A3)-(fin*(q-p))), ipx7);
		     z=0;
		     for(ii=fin*(q-p);ii<length(A3);ii++){
			     REAL(auxA3)[z]=REAL(A3)[ii];
			     z=z+1;
		     }
		     REPROTECT(A3 = allocVector(REALSXP,  length(auxA3)), ipx6);
		     for(ii=0;ii<length(auxA3);ii++)REAL(A3)[ii]=REAL(auxA3)[ii];
	    }	 
    } 
  /*  SEXP m3;
    ncol=q-p; 
    nrow=length(A3)/(q-p);
    PROTECT(m3 = allocMatrix(REALSXP, nrow, ncol));
    ++nProtected;
    for(ii=0;ii<nrow;ii++){
	    for(jj=0;jj<ncol;jj++){
		    REAL(m3)[ii+(nrow*jj)]=REAL(A3)[ii*ncol+jj];
	    }
    } */

    SEXP mrl3;
    PROTECT(mrl3 = allocVector(REALSXP, q-p));
    ++nProtected;
    z=0;
    for(ii=p; ii<q;ii++){
	    REAL(mrl3)[z]=REAL(mrl)[ii];
	    z++;
    }
    SEXP means3;
    PROTECT(means3 = allocMatrix(REALSXP, q-p, q-p));
    ++nProtected;
    y=0;
    for(ii=p;ii<q;ii++){
	    z=0;
	    for(jj=p;jj<q;jj++){
		    REAL(means3)[y+(q-p)*z]=REAL(means)[ii+length(data)*jj];
		    z++;
	    }
	    y++;
    } 
/*	      Rprintf("\n\nnA3=");
	     for(ii=0;ii<length(A3);ii++)Rprintf(" %f, ",REAL(A3)[ii]);

	     Rprintf("\nmrl3 length= %d \n", length(mrl3)); */

	     solution=cirPAVA(G3, A3, means3, mrl3, 3);
	     
	  /*   Rprintf("\nsolution=");
	     for(ii=0;ii<length(solution);ii++)Rprintf(" %f, ",REAL(solution)[ii]); */

    SET_VECTOR_ELT(lA3, p, solution);
   } 


/* RESULTS: */

   SET_VECTOR_ELT(result, 0, lA1);
   SET_VECTOR_ELT(result, 1, lA2);
   SET_VECTOR_ELT(result, 2, lA3);

   UNPROTECT(nProtected);
   
   return result;
}


			


