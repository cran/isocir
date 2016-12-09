#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP cirPAVA(SEXP, SEXP, SEXP, SEXP, int);


SEXP funG1(SEXP data, int s, SEXP means, SEXP mrl){
	
	int i, j, jj, exit, a, b, z, ii, w, z2, nProtected = 0, fin, d, f, count, nr, add;
	double soli, check, aux=0.0;
   	
	SEXP A1;
	SEXP G1;
	SEXP solution;
	SEXP auxA1;
	
	/* G1 */
     	
	
     	PROTECT_INDEX ipx10;
     	PROTECT_WITH_INDEX(G1 = allocVector(REALSXP,  s+1), &ipx10);
     	++nProtected;
      	
     	PROTECT_INDEX ipx4;
     	PROTECT_WITH_INDEX(A1 = allocVector(REALSXP,  s+1), &ipx4);
     	++nProtected;
	
     	PROTECT_INDEX ipx2;
     	PROTECT_WITH_INDEX(solution = allocVector(REALSXP,  s+1), &ipx2);
     	++nProtected;
	
     	PROTECT_INDEX ipx7;
     	PROTECT_WITH_INDEX(auxA1 = allocVector(REALSXP,  s+1), &ipx7);
     	++nProtected;
	 
	
     	
    	for(ii=0;ii<length(A1);ii++)REAL(A1)[ii]=REAL(data)[ii];
    	for(ii=0;ii<length(G1);ii++)REAL(G1)[ii]=REAL(data)[ii];
	
for(i=s;i>=0;i--){
		   
       		   fin=(length(A1)/(s+1));
       		   for(j=0;j<fin;j++){
       			   aux=REAL(A1)[j*(s+1)+i];
	  		   
       			   if(aux<=M_PI){
       				   z=length(A1);
       				   REPROTECT(auxA1 = allocVector(REALSXP,  length(A1)), ipx7);
       				   for(ii=0;ii<length(A1);ii++)REAL(auxA1)[ii]=REAL(A1)[ii];
       				   
       				   REPROTECT(A1 = allocVector(REALSXP,  z+(s+1)), ipx4);
       				   for(ii=0;ii<length(auxA1);ii++)REAL(A1)[ii]=REAL(auxA1)[ii];
       				   for(ii=0;ii<s+1;ii++){
       					   REAL(A1)[z]=REAL(auxA1)[j*(s+1)+ii];
       					   z=z+1;
       				   }
       			   }		
       			   if(aux<=(2*M_PI) && aux>=(3*M_PI/2)){
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
   				   /*  && REAL(solution)[0]!=0 */
       				   if(add==1){
       					   z=length(A1);
       					   REPROTECT(auxA1 = allocVector(REALSXP,  length(A1)), ipx7);
       					   for(ii=0;ii<length(A1);ii++)REAL(auxA1)[ii]=REAL(A1)[ii];
       					   
       					   REPROTECT(A1 = allocVector(REALSXP,  z+(s+1)), ipx4);
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
       					   REPROTECT(auxA1 = allocVector(REALSXP,  length(A1)), ipx7);	   
       					   for(ii=0;ii<length(A1);ii++)REAL(auxA1)[ii]=REAL(A1)[ii];
       					   REPROTECT(A1 = allocVector(REALSXP,  z+(s+1)), ipx4);
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
      						   REPROTECT(auxA1 = allocVector(REALSXP,  length(A1)), ipx7);
      						   for(ii=0;ii<length(A1);ii++)REAL(auxA1)[ii]=REAL(A1)[ii];
      						   
      						   REPROTECT(A1 = allocVector(REALSXP,  z+(s+1)), ipx4);
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
      							   REPROTECT(auxA1 = allocVector(REALSXP,  length(A1)), ipx7);
      							   for(ii=0;ii<length(A1);ii++)REAL(auxA1)[ii]=REAL(A1)[ii];
      							   REPROTECT(A1 = allocVector(REALSXP,  z+(s+1)), ipx4);
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
      						   if(add==1 ){
      							   z=length(A1);
      							   REPROTECT(auxA1 = allocVector(REALSXP,  length(A1)), ipx7);
      							   for(ii=0;ii<length(A1);ii++)REAL(auxA1)[ii]=REAL(A1)[ii];
      							   REPROTECT(A1 = allocVector(REALSXP,  z+(s+1)), ipx4);
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
       			   w=0;
       			   for(ii=0;ii < ((length(A1)/(s+1))-fin); ii++){
       				   if(REAL(A1)[(fin+ii)*(s+1)]!=0)w=w+1;
       			   }
       			   if(w > 0){
       				   REPROTECT(auxA1 = allocVector(REALSXP, w*(s+1)), ipx7);
       				   z=0;
       				   for(ii=0;ii < ((length(A1)/(s+1))-fin); ii++){
       					   if(REAL(A1)[(fin+ii)*(s+1)]!=0){
       						   for(z2=0;z2<(s+1);z2++){
       							   REAL(auxA1)[z]=REAL(A1)[(fin+ii)*(s+1)+z2];
       							   z=z+1;
       						   }
       					   }
       				   }
       			   }
       			   if(w==0){
      				   REPROTECT(auxA1 = allocVector(REALSXP,(s+1)), ipx7);
      				   for(ii=0;ii<(s+1);ii++){				   
					   REAL(auxA1)[ii]=9999;
					   
      				   }
       			   }
       			   REPROTECT(A1 = allocVector(REALSXP,  length(auxA1)), ipx4);
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

	 if(REAL(A1)[0]>7){
       		   for(ii=0;ii<(s+1);ii++)REAL(solution)[ii]=REAL(A1)[ii];
       	   }

    	if(REAL(A1)[0]<7){
		
	    	SEXP means1;
	    	PROTECT(means1 = allocMatrix(REALSXP, s+1, s+1));
	    	++nProtected;
    		for(ii=0;ii<s+1;ii++){
	    		for(jj=0;jj<s+1;jj++)REAL(means1)[ii+(s+1)*jj]=REAL(means)[ii+length(data)*jj];
	    	} 
	     	solution=cirPAVA(G1, A1, means1, mrl, 1);
	}
	UNPROTECT(nProtected);
	return(solution);
}

SEXP funG2(SEXP data, int p, int s, SEXP means, SEXP mrl){

		
	R_len_t q=length(data);

	int i, j, jj, y, L, exit,a, b, z, ii, nProtected = 0, fin, d, f, count, nr, add;
        double soli, check, aux=0.0;

	L=p-s;

	SEXP A2;
	SEXP G2;
	SEXP solution;
	SEXP auxA2;

/* G2 */

   PROTECT_INDEX ipx2;
   PROTECT_WITH_INDEX(solution = allocVector(REALSXP,  L), &ipx2);
   ++nProtected;

   PROTECT_INDEX ipx5;
   PROTECT_WITH_INDEX(A2 = allocVector(REALSXP, L), &ipx5);
   ++nProtected;

   PROTECT_INDEX ipx8;
   PROTECT_WITH_INDEX(auxA2 = allocVector(REALSXP,  L), &ipx8);
   ++nProtected;

    PROTECT_INDEX ipx11;
   PROTECT_WITH_INDEX(G2 = allocVector(REALSXP,  L), &ipx11);
   ++nProtected;

     

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
	     for(i=0;i<L;i++){		  
		     fin=length(A2)/L;
		     for(j=0;j<fin;j++){			    
			     aux=REAL(A2)[j*L+i];			    		     
			     if(aux>=M_PI/2 && aux<=3*M_PI/2){
				     for(ii=0;ii<L;ii++)REAL(solution)[ii]=REAL(A2)[j*L+ii];

				     z=length(A2);
				     REPROTECT(auxA2 = allocVector(REALSXP,  length(A2)), ipx8);
				     for(ii=0;ii<length(A2);ii++)REAL(auxA2)[ii]=REAL(A2)[ii];
				     REPROTECT(A2 = allocVector(REALSXP,  z+L), ipx5);
				     for(ii=0;ii<length(auxA2);ii++)REAL(A2)[ii]=REAL(auxA2)[ii];
				     for(ii=0;ii<length(solution);ii++){
					     REAL(A2)[z]=REAL(solution)[ii];
					     z=z+1;
				     }

			     }
			     if(aux>=3*M_PI/2 || aux<=M_PI/2){
				     b=i;
				     a=i;
				     z=L-1;
				     check=REAL(A2)[(j*L)+a];
				     while(check!=REAL(A2)[(j*L)+z])z--;
				     a=z;
				     
				     soli=aux;
				     exit=0;
				    
				     while((soli>3*M_PI/2 || soli<M_PI/2) && exit==0){
					     b=b-1;
					     if(b<0)exit=1;
					     if(b>=0){
						     z=0;
						     check=REAL(A2)[(j*L)+b];
						     while(check!=REAL(A2)[(j*L)+z])z++;
						     b=z;
						     
						     soli=REAL(means)[(b+s)+length(data)*(a+s)];
					     }
				     }
				     
				     if(exit==0){
					     for(ii=0;ii<(p-s);ii++)REAL(solution)[ii]=REAL(A2)[j*L+ii];
					     for(ii=b;ii<=a;ii++)REAL(solution)[ii]=soli;

					      add=1;
		  			      nr=(length(A2)-(fin*L))/L;
		  			      if(length(A2)>fin*L){
		  				      for(d=0;d<nr;d++){
		  					      count=0;
		  					      for(f=0;f<L;f++){
		  						      if(REAL(A2)[((d+fin)*L)+f]==REAL(solution)[f])count++;
		  					      }
		  					      if(count==L){
		   						      add=0;
		   						      break;
		  					      }
		  				      }
		  			      }
		  			      if(add==1){
	 					      z=length(A2);
	 					      REPROTECT(auxA2 = allocVector(REALSXP,  length(A2)), ipx8);
	 					      for(ii=0;ii<length(A2);ii++)REAL(auxA2)[ii]=REAL(A2)[ii];
	 					      REPROTECT(A2 = allocVector(REALSXP,  z+L), ipx5);
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
				     check=REAL(A2)[(j*L)+a];
				     while(check!=REAL(A2)[(j*L)+z])z++;
				     a=z;
				     soli=aux;
				     exit=0;
				     while((soli>3*M_PI/2 || soli<M_PI/2) && exit==0){
					     b=b+1;
					     if(b>=L)exit=1;
					     if(b<L){
						     z=L-1;
						     check=REAL(A2)[(j*L)+b];
						     while(check!=REAL(A2)[(j*L)+z])z--;
						     b=z;
						     
						     soli=REAL(means)[(a+s)+length(data)*(b+s)];
					     }
				     }
				     
				     if(b==L)b=L-1;
				     if(exit==0){
					     for(ii=0;ii<L;ii++)REAL(solution)[ii]=REAL(A2)[j*L+ii];
					     for(ii=a;ii<=b;ii++)REAL(solution)[ii]=soli;

					      add=1;
					      nr=(length(A2)-(fin*L))/L;
		  			      if(length(A2)>fin*L){
		  				      for(d=0;d<nr;d++){
		  					      count=0;
		  					      for(f=0;f<L;f++){
		  						      if(REAL(A2)[((d+fin)*L)+f]==REAL(solution)[f])count++;
		  					      }
		  					      if(count==L){
		   						      add=0;
		   						      break;
		  					      }
		  				      }
		  			      }
		  			      if(add==1){					     
					     z=length(A2);
					     REPROTECT(auxA2 = allocVector(REALSXP,  length(A2)), ipx8);
					     for(ii=0;ii<length(A2);ii++)REAL(auxA2)[ii]=REAL(A2)[ii];
					     
					     REPROTECT(A2 = allocVector(REALSXP,  z+L), ipx5);
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
					     check=REAL(A2)[(j*L)+a];
					     while(check!=REAL(A2)[(j*L)+z])z++;
					     a=z;
					     
					     b=i+1;
					     z=L-1;
					     check=REAL(A2)[(j*L)+b];
					     while(check!=REAL(A2)[(j*L)+z])z--;
					     b=z;
					    
					     for(ii=0;ii<L;ii++)REAL(solution)[ii]=REAL(A2)[j*L+ii];
					     for(ii=a;ii<=b;ii++)REAL(solution)[ii]=REAL(means)[(a+s)+length(data)*(b+s)];	     
					    
					      add=1;
		  			      nr=(length(A2)-(fin*L))/L;
		  			      if(length(A2)>fin*L){
		  				      for(d=0;d<nr;d++){
		  					      count=0;
		  					      for(f=0;f<L;f++){
		  						      if(REAL(A2)[((d+fin)*L)+f]==REAL(solution)[f])count++;
		  					      }
		  					      if(count==L){
		   						      add=0;
		   						      break;
		  					      }
		  				      }
		  			      }
		  			      if(add==1){				     				     
					     z=length(A2);
					     REPROTECT(auxA2 = allocVector(REALSXP,  length(A2)), ipx8);
					     for(ii=0;ii<length(A2);ii++)REAL(auxA2)[ii]=REAL(A2)[ii];
					     REPROTECT(A2 = allocVector(REALSXP,  z+L), ipx5);
					     for(ii=0;ii<length(auxA2);ii++)REAL(A2)[ii]=REAL(auxA2)[ii];
					     for(ii=0;ii<length(solution);ii++){
						     REAL(A2)[z]=REAL(solution)[ii];
						     z=z+1;
					     } 
					      }
				     }

			     } 
		     }
		    
		     if((length(A2)/L)>fin){
			     REPROTECT(auxA2 = allocVector(REALSXP,  length(A2)-(fin*L)), ipx8);
			     z=0;
			     for(ii=fin*L;ii<length(A2);ii++){
				     REAL(auxA2)[z]=REAL(A2)[ii];
				     z=z+1;
			     }
			     REPROTECT(A2 = allocVector(REALSXP,  length(auxA2)), ipx5);
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
	  /*   Rprintf("\n antes");*/
/*	     for(ii=0;ii<length(solution);ii++)Rprintf(" %f, ",REAL(A2)[ii]); */
	     solution=cirPAVA(G2, A2, means2, mrl2, 2);
	     UNPROTECT(nProtected);
/*Rprintf("\n despues");
	      for(ii=0;ii<length(solution);ii++)Rprintf(" %f, ",REAL(solution)[ii]);*/
	     return(solution);
}


SEXP funCIREnuevo(SEXP data, SEXP means, SEXP mrl){
	
	R_len_t q=length(data);
     	int i, j, jj, y, exit, s, p, a, b, z, ii, w, z2, nProtected = 0, fin, d, f, count, nr, add, ncol, nrow, condition;
     	double soli, check, aux=0.0, bsce, sce0;
     
	SEXP lA3;
     	SEXP solution;
   	SEXP csolution;
     	SEXP Msolution;

 /*  Rprintf("data: ");
    for(ii=0;ii<q;ii++)Rprintf(" %f, ",REAL(data)[ii]);
    Rprintf("\n");*/
 

	 /* G3 */
	
	PROTECT(lA3 = allocVector(VECSXP, q));
     	++nProtected;
	
     	SEXP A3, auxA3, G3;
	
     	PROTECT_INDEX ipx6;
     	PROTECT_WITH_INDEX(A3 = allocVector(REALSXP,  0), &ipx6);
     	++nProtected;
	
     	PROTECT_INDEX ipx9;
     	PROTECT_WITH_INDEX(auxA3 = allocVector(REALSXP,  0), &ipx9);
     	++nProtected;
     	
     	PROTECT_INDEX ipx12;
     	PROTECT_WITH_INDEX(G3 = allocVector(REALSXP,  0), &ipx12);
     	++nProtected;
	
    	PROTECT_INDEX ipx2;
     	PROTECT_WITH_INDEX(solution = allocVector(REALSXP, 0), &ipx2);
     	++nProtected;
	

	
	
	
     	for(p=0;p<q;p++){
		
		REPROTECT(A3 = allocVector(REALSXP, q-p), ipx6);
	    	i=0;
	    	for(ii=p;ii<q;ii++){
		   	REAL(A3)[i]=REAL(data)[ii];
		   	i=i+1;
	    	}
	    	REPROTECT(G3 = allocVector(REALSXP, q-p), ipx12);
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
		    if(aux>=M_PI){
			    z=length(A3);
		
			    REPROTECT(auxA3 = allocVector(REALSXP,  length(A3)), ipx9);
			    for(ii=0;ii<length(A3);ii++)REAL(auxA3)[ii]=REAL(A3)[ii];

                            REPROTECT(A3 = allocVector(REALSXP,  z+(q-p)), ipx6);
			    for(ii=0;ii<length(auxA3);ii++)REAL(A3)[ii]=REAL(auxA3)[ii];
			    for(ii=0;ii<(q-p);ii++){
				    REAL(A3)[z]=REAL(auxA3)[j*(q-p)+ii];
				    z=z+1;
			    }
		    }
		    if(aux<=M_PI/2){
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
			 
			    REPROTECT(auxA3 = allocVector(REALSXP,  length(A3)), ipx9);
			    for(ii=0;ii<length(A3);ii++)REAL(auxA3)[ii]=REAL(A3)[ii];

                            REPROTECT(A3 = allocVector(REALSXP,  z+(q-p)), ipx6);

			    for(ii=0;ii<length(auxA3);ii++)REAL(A3)[ii]=REAL(auxA3)[ii];
			    for(ii=0;ii<q-p;ii++){
				    REAL(A3)[z]=REAL(solution)[ii];
				    z=z+1;
			    }
			    }
		    }
		    if(aux>=M_PI/2 && aux<=M_PI){

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
				  
				    REPROTECT(auxA3 = allocVector(REALSXP,  length(A3)), ipx9);
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
				    REPROTECT(auxA3 = allocVector(REALSXP,  z), ipx9);
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
					    REPROTECT(auxA3 = allocVector(REALSXP,  length(A3)), ipx9);
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
					    REPROTECT(auxA3 = allocVector(REALSXP,  length(A3)), ipx9);
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
		    w=0;
		    for(ii=0;ii < ((length(A3)/(q-p))-fin); ii++){
			    if(REAL(A3)[((fin+ii)*(q-p))+(q-p-1)]!=2*M_PI)w=w+1;
		    }
		    if(w > 0){
			    REPROTECT(auxA3 = allocVector(REALSXP, w*(q-p)), ipx9);
			    z=0;
			    for(ii=0;ii < ((length(A3)/(q-p))-fin); ii++){
				    if(REAL(A3)[(fin+ii)*(q-p)+(q-p-1)]!=2*M_PI){
					    for(z2=0;z2<(q-p);z2++){
						    REAL(auxA3)[z]=REAL(A3)[(fin+ii)*(q-p)+z2];
						    z=z+1;
					    }
				    }
			    }
		    }
		    if(w==0){
			     REPROTECT(auxA3 = allocVector(REALSXP,(q-p)), ipx9);
			     for(ii=0;ii<(q-p);ii++){
				     REAL(auxA3)[ii]=9999;
			     }
		    }
		    REPROTECT(A3 = allocVector(REALSXP,  length(auxA3)), ipx6);
		    for(ii=0;ii<length(auxA3);ii++)REAL(A3)[ii]=REAL(auxA3)[ii];

	    }	 
    } 
  /*  SEXP m3;
     int nrow, ncol;
    ncol=q-p; 
    nrow=length(A3)/(q-p);
    PROTECT(m3 = allocMatrix(REALSXP, nrow, ncol));
    ++nProtected;
    for(ii=0;ii<nrow;ii++){
	    for(jj=0;jj<ncol;jj++){
		    REAL(m3)[ii+(nrow*jj)]=REAL(A3)[ii*ncol+jj];
	    }
    } */
	
      		if(REAL(A3)[0]>7){
	    		for(ii=0;ii<(q-p);ii++)REAL(solution)[ii]=REAL(A3)[ii];
	    	}
	    	if(REAL(A3)[0]<7){
			
			
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
		   	solution=cirPAVA(G3, A3, means3, mrl3, 3);
	    	}
	    	SET_VECTOR_ELT(lA3, p, solution);
		/*    SET_VECTOR_ELT(lA3, p, m3);*/
     	} 
	 /*	Rprintf("fin G3");*/

/* Aqu? empieza la union de los tres cachos */  
	
       add=0;
       
       SEXP A1;
       SEXP A2;
       SEXP auxM;
       SEXP porsi;
       
       PROTECT_INDEX ipx13;
       PROTECT_WITH_INDEX(auxM = allocVector(REALSXP, 0), &ipx13);
       ++nProtected;
       
       PROTECT_INDEX ipx3;
       PROTECT_WITH_INDEX(Msolution = allocVector(REALSXP, 0), &ipx3);
       ++nProtected;
       
       PROTECT_INDEX ipx4;
       PROTECT_WITH_INDEX(A1 = allocVector(REALSXP, 0), &ipx4);
       ++nProtected;
        PROTECT_INDEX ipx15;
       PROTECT_WITH_INDEX(porsi = allocVector(REALSXP, 0), &ipx15);
       ++nProtected;
       
       PROTECT_INDEX ipx5;
       PROTECT_WITH_INDEX(A2 = allocVector(REALSXP, 0), &ipx5);
       ++nProtected;

       PROTECT_INDEX ipx1;
     	PROTECT_WITH_INDEX(csolution = allocVector(REALSXP, 0), &ipx1);
     	++nProtected;

/*	REPROTECT(A1 = allocVector(REALSXP, 1), ipx4);
   	REAL(A1)[0]=999.9;*/
 /*	Rprintf("entrando a bucle s con if");*/
       for(s=-1;s<q;s++){
	     /*  Rprintf("s= %d, ",s);*/
	       if(s> -1){
   		       REPROTECT(A1 = allocVector(REALSXP, s+1), ipx4);
		        REPROTECT(porsi = allocVector(REALSXP, s+1), ipx15);
   		       A1=funG1(data, s, means, mrl);
		     /*  Rprintf("A1: ");
		       for(ii=0;ii<length(A1);ii++)Rprintf(" %f, ",REAL(A1)[ii]);					
		       Rprintf("hecho G1");*/
		       for(ii=0;ii<length(A1);ii++)REAL(porsi)[ii]=REAL(A1)[ii];					
      		    /* for(ii=0;ii<length(porsi);ii++)Rprintf(" %f, ",REAL(porsi)[ii]);*/
   	       }
    	       if(s==-1 || REAL(A1)[0]<7){
   		       for(p=s+1;p<(q+1);p++){
			      /* Rprintf("p= %d, ",p);*/
    			       condition=0;
			       if(p<q){
				       if(REAL(VECTOR_ELT(lA3, p))[0]<7)condition=1;
			       }
			       if(p==q)condition=1;
			       if(condition==1){
				       REPROTECT(csolution = allocVector(REALSXP, q), ipx1);
				      /*  Rprintf("entrando");*/
				       if(p<q){
					       REPROTECT(A3 = allocVector(REALSXP, q-p), ipx6);
					       for(ii=0;ii<length(A3);ii++)REAL(A3)[ii]=REAL(VECTOR_ELT(lA3, p))[ii];
				       }
				       if(p==(s+1)){
					       if(s==-1){						       
						       for(ii=0;ii<q;ii++)REAL(csolution)[ii]=REAL(A3)[ii];
						       add=1;
					       }
					       
					       if(s>-1 && s<(q-1)){						      
   						       for(ii=0;ii<length(csolution);ii++){
   							       if(ii<s+1)REAL(csolution)[ii]=REAL(porsi)[ii];
   							       if(ii>=p)REAL(csolution)[ii]=REAL(A3)[ii-p];
   						       }
						       add=1;
   					       }
   					       if(s==(q-1)){						      
  						       for(ii=0;ii<q;ii++)REAL(csolution)[ii]=REAL(porsi)[ii];
  						       add=1;
   					       }
   				       }

   				       if(p>(s+1)){
					       REPROTECT(A2 = allocVector(REALSXP, p-s-1), ipx5);
					       A2=funG2(data, p, s+1, means, mrl);
					    /*   Rprintf("A2: ");
					       for(ii=0;ii<length(A2);ii++)Rprintf(" %f, ",REAL(A2)[ii]);
					       Rprintf("hecho G2");*/
					       if(s==-1 && p<q){						      
						       	for(ii=0;ii<length(csolution);ii++){
 								if(ii<p)REAL(csolution)[ii]=REAL(A2)[ii];
 								if(ii>=p)REAL(csolution)[ii]=REAL(A3)[ii-p];
 							}
 						 /*   add=1;*/
 						}					      
 						if(s==-1 && p==q){       						
 							for(ii=0;ii<q;ii++)REAL(csolution)[ii]=REAL(A2)[ii];
                                                     /* Rprintf("que se convierte en: %f", REAL(csolution)[0]);*/
 						 /*   add=1;*/
 						}
 						if(s>-1 && p<q){
						/*	Rprintf("Comprobar A1 con porsi: ");
							for(ii=0;ii<length(porsi);ii++)Rprintf(" %f, ",REAL(porsi)[ii]);					
							Rprintf("Comprobar A1: ");
	 						for(ii=0;ii<length(A1);ii++)Rprintf(" %f, ",REAL(A1)[ii]);					       */

						       for(ii=0;ii<length(csolution);ii++){
							       if(ii<s+1)REAL(csolution)[ii]=REAL(porsi)[ii];
   							       if(ii>s && ii<p)REAL(csolution)[ii]=REAL(A2)[ii-(s+1)];
   							       if(ii>=p)REAL(csolution)[ii]=REAL(A3)[ii-p];
   						       }
						        /*   add=1;*/
 						}
 						if(s>-1 && p==q){
							/*		Rprintf("Comprobar A1 con porsi: ");
							for(ii=0;ii<length(porsi);ii++)Rprintf(" %f, ",REAL(porsi)[ii]);
						Rprintf("Comprobar A1: ");
	 						for(ii=0;ii<length(A1);ii++)Rprintf(" %f, ",REAL(A1)[ii]);					       
								Rprintf("entra aqui");*/
							for(ii=0;ii<length(csolution);ii++){
 								if(ii<(s+1))REAL(csolution)[ii]=REAL(porsi)[ii];
 								if(ii>s && ii<p)REAL(csolution)[ii]=REAL(A2)[ii-(s+1)];
							}
								
 						 /*   add=1;*/
						}
						aux=0;
						for(ii=0;ii<length(csolution);ii++)aux=aux+REAL(csolution)[ii];
 						if(aux<1000)add=1;
            				/*	Rprintf("aux= %f", aux);*/
				       }				       					       
			       }
			       if(add==1){
				       add=0;
				       z=length(Msolution);
				       if(z>0){
				       REPROTECT(auxM = allocVector(REALSXP,  z), ipx13);
				       for(ii=0;ii<z;ii++)REAL(auxM)[ii]=REAL(Msolution)[ii];
				       }
				       REPROTECT(Msolution = allocVector(REALSXP,  z+q), ipx3);
				       if(z>0){
					       for(ii=0;ii<z;ii++)REAL(Msolution)[ii]=REAL(auxM)[ii];
				       }
				       for(ii=z;ii<z+q;ii++)REAL(Msolution)[ii]=REAL(csolution)[ii-z];
				       
			       	       REPROTECT(csolution = allocVector(REALSXP, 0), ipx1);
			       }
		       }
   	       }   	      
       }
       if(length(Msolution)==0){
	        REPROTECT(Msolution = allocVector(REALSXP,  q), ipx3);
		for(ii=0;ii<q;ii++)REAL(Msolution)[ii]=9999;
       }
      /* Rprintf("Msolution: ");
       for(ii=0;ii<length(Msolution);ii++)Rprintf(" %f, ",REAL(Msolution)[ii]);*/

    /*Convertir Msolution en matriz*/
    /*    SEXP salida; 
       ncol=q;
       nrow=length(Msolution)/q;
       PROTECT(salida = allocMatrix(REALSXP, nrow, ncol));
       ++nProtected;
       for(ii=0;ii<nrow;ii++){
	       for(jj=0;jj<ncol;jj++)REAL(salida)[ii+(nrow*jj)]=REAL(Msolution)[ii*ncol+jj];
       }*/
    

       SEXP CIRE; 
       SEXP aux2;
       PROTECT(CIRE =  allocVector(REALSXP, q+1));
       ++nProtected;
       PROTECT(aux2 =  allocVector(REALSXP, q));
       ++nProtected;
       ncol=q;
       nrow=length(Msolution)/q;
       bsce=9999;
       for(ii=0;ii<nrow;ii++){
	       for(jj=0;jj<q;jj++)REAL(aux2)[jj]=REAL(Msolution)[ii*ncol+jj];
	       sce0=0;
	       	for(y=0;y<length(data);y++)sce0=sce0+(REAL(mrl)[y]*(1-cos(REAL(data)[y]-REAL(aux2)[y])));
	      /*  Rprintf("sce= %f, ",sce0);*/
	       if(sce0<bsce){
		       for(z2=0;z2<=q;z2++){
			       if(z2<q)REAL(CIRE)[z2]=REAL(aux2)[z2];
			       if(z2==q)REAL(CIRE)[z2]=sce0;
		       }
		       bsce=sce0;
	       }
       }
       
    
       UNPROTECT(nProtected);

     /*  Rprintf(" sale del todo. \n");*/
    
       return (CIRE);
}

