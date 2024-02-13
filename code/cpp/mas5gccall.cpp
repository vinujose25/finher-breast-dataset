#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <typeinfo>
#include <cstddef>
#include <string>
#include <vector>
#include <string>

#include "./alglib-3.8.2/src/statistics.h" 
#include "dabg2call.h"

//.Call("mas5gccall",data,dataNcMed,indexList,gcList)

extern "C" SEXP mas5gccall(SEXP intensity, SEXP ncMedIntensity, 
						       SEXP indexList, SEXP gcList, SEXP tau) {
	
	ExprMat expr(intensity); //REprintf("Done sampMat.\n");
	ExprMat ncMed(ncMedIntensity); //REprintf("Done sampMat.\n");
	ProbeSets psets(indexList,gcList); //REprintf("Done probesets.\n");
	
	double tau2=*REAL(tau); 
	std::vector<int> prIndex, gcIndex;
	int prIndexSize, prSetCount=psets.get_prSetCount();
	alglib::real_1d_array pr;
	double btail, ltail, rtail; //return pvalues
	SEXP dCal=PROTECT(allocMatrix(REALSXP,prSetCount,expr.get_nCol())); // for storing detection pval
	std::vector<double> curSamp, curNcMed; //current sample
	curSamp.reserve(expr.get_nRow());
        curNcMed.reserve(ncMed.get_nRow());
	
	for(int i=0;i<expr.get_nCol();++i){
		REprintf(expr.get_colName(i)); PrintValue(ScalarInteger(i)); 
		curSamp.clear();
		curNcMed.clear();
		expr.get_column(i,curSamp);	// return by reference; &curSamp
		ncMed.get_column(i,curNcMed);  // return by reference; &curNcMed
		for(int j=0;j<prSetCount;++j){
//			if(j%1000==0){PrintValue(ScalarInteger(j));}
			psets.get_prIndex(j,prIndex);
			psets.get_refIndex(j,gcIndex);
			prIndexSize=prIndex.size();
			
			pr.setlength(prIndexSize); 
			for(int k=0;k<prIndexSize;++k){ // pr
				pr[k]=curSamp[prIndex[k]-1] - curNcMed[gcIndex[k]-3];
			}

			if(prIndexSize>=5){
				alglib::wilcoxonsignedranktest(pr,prIndexSize, tau2, btail,ltail,rtail);
			}else{
				rtail=99.99; 	
			}
			REAL(dCal)[ (i*prSetCount)+j ] = rtail;  	
		} // end probesets loop							
	} // end expr loop
	
	UNPROTECT(1); // unprotect dCal		
	return(dCal);

} // end extern C		


