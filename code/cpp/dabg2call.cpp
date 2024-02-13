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

//.Call("dabg2call",data,indexList,indexList_ref2,refProbeSetType)

extern "C" SEXP dabg2call(SEXP intensity, SEXP indexList, SEXP indexList_ref,
						 SEXP refProbeSetType) {
	
	std::string bgType(CHAR(asChar(refProbeSetType)));
	if(bgType=="medGC" || bgType=="maxGC" || bgType=="relGC"){
//		REprintf((bgType+"\n").data());
		if(length(indexList)!=length(indexList_ref)){
			REprintf("Non congruant gc and index lists\n");
		}
	}else if(bgType=="allGC" || bgType=="fixGC"){
//		REprintf((bgType+"\n").data());
		if(length(indexList)==length(indexList_ref)){
			REprintf("More than one element in indexList_ref\n");
		}		
	}
	
	ExprMat expmat(intensity); //REprintf("Done sampMat.\n");
	ProbeSets psets(indexList,indexList_ref); //REprintf("Done probesets.\n");
	
	std::vector<int> prIndex, refIndex;
	int prIndexSize,refIndexSize, prSetCount=psets.get_prSetCount();
	alglib::real_1d_array pr, ref;
	double btail, ltail, rtail; //return varibales
	SEXP dCal=PROTECT(allocMatrix(REALSXP,prSetCount,expmat.get_nCol())); // for detection pval
	std::vector<double> curSamp; //current sample
	curSamp.reserve(expmat.get_nRow()); 
	
	if(bgType=="medGC" || bgType=="maxGC" || bgType=="relGC"){
		
//		REprintf("medGc|maxGc|relGc\n");
		for(int i=0;i<expmat.get_nCol();++i){
			REprintf(expmat.get_colName(i));  PrintValue(ScalarInteger(i)); 
			curSamp.clear();
			expmat.get_column(i,curSamp);	// return by reference; &curSamp
			for(int j=0;j<prSetCount;++j){
//				if(j%1000==0){PrintValue(ScalarInteger(j));}
				psets.get_prIndex(j,prIndex);
				psets.get_refIndex(j,refIndex);
				prIndexSize=prIndex.size();
				refIndexSize=refIndex.size();
				
				pr.setlength(prIndexSize); //wilcox test function
				ref.setlength(refIndexSize);
				for(int k=0;k<prIndexSize;++k){ // pr
					pr[k]=curSamp[prIndex[k]-1];
				}
				for(int k=0;k<refIndexSize;++k){ // ref
					ref[k]=curSamp[refIndex[k]-1];
				}
	
				if(prIndexSize>=5){
					alglib::mannwhitneyutest(pr,prIndexSize,ref,refIndexSize,btail,ltail,rtail);					
				}else{
					rtail=99.99; 	
				}
				REAL(dCal)[ (i*prSetCount)+j ] = rtail;  	
			} // end probesets loop							
		} // end expmat loop

	}else if(bgType=="allGC" || bgType=="fixGC"){

//		REprintf("allGc|fixGc\n");
		psets.get_refIndex(0,refIndex);
		refIndexSize=refIndex.size();
		ref.setlength(refIndexSize);
		
		for(int i=0;i<expmat.get_nCol();++i){
			REprintf(expmat.get_colName(i));  PrintValue(ScalarInteger(i)); 
			curSamp.clear();
			expmat.get_column(i,curSamp);	// return by reference; &curSamp
			
			for(int k=0;k<refIndexSize;++k){ // ref for each sample
				ref[k]=curSamp[refIndex[k]-1];
			}

			for(int j=0;j<prSetCount;++j){ //pr for each probeSet
//				if(j%1000==0){PrintValue(ScalarInteger(j));}
				psets.get_prIndex(j,prIndex);
				prIndexSize=prIndex.size();
				pr.setlength(prIndexSize); //wilcox test function
				for(int k=0;k<prIndexSize;++k){ // pr
					pr[k]=curSamp[prIndex[k]-1];
				}
				
				if(prIndexSize>=5){
					alglib::mannwhitneyutest(pr,prIndexSize,ref,refIndexSize,btail,ltail,rtail); 	
				}else{
					rtail=99.99;
				}
				REAL(dCal)[ (i*prSetCount)+j ] = rtail;  	
			} // end probesets loop							
		} // end expmat loop

	}
	
	UNPROTECT(1); // unprotect dCal		
	return(dCal);

} // end extern C		


