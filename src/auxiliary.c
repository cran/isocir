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
		A=A/length(data);
		B=B/length(data);
	}
	if ((A >= 0) && (B > 0))result=atan(A/B);
        if ((A > 0) && (B == 0))result=(M_PI)/2;
    	if (B < 0)result=atan(A/B) + M_PI;
	if ((A < 0) && (B >= 0))result=atan(A/B) + (2 * M_PI);
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

	PROTECT_INDEX ipx14;
	PROTECT_WITH_INDEX(aux3 = allocVector(REALSXP,  0), &ipx14);
     	++nProtected;
        

	int nrow, ii, i, j, sumproblems, z, a, b;
	double check, Sumce=9999, nSumce, condition;
	
	nrow=length(A)/len;
	for(j=0;j<len;j++)REAL(result)[j]=9999;

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
	
		REPROTECT(aux3 = allocVector(REALSXP,  sumproblems), ipx14);
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
			REPROTECT(aux3 = allocVector(REALSXP,  sumproblems), ipx14);
			for(j=0;j<length(aux3);j++)REAL(aux3)[j]=REAL(aux2)[j];
		} 
		nSumce=sce(aux, datas, mrl);
	/*	Rprintf("SCE: %f",nSumce);*/
		if(group==1)condition5=REAL(aux)[length(aux)-1]<=(M_PI/2);
		if(group==2)condition5=(REAL(aux)[0]>M_PI/2)&&(REAL(aux)[length(aux)-1]<=(3*M_PI/2));
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

