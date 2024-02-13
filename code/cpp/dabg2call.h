// getDABG2.h

class ExprMat{
	private:	
	SEXP data;
	int * dim; // change it to ncol and nrow
	std::vector<const char *> colnames;
	std::vector<const char *> rownames;
	public:	
	ExprMat(SEXP);
	int get_nRow();
	int get_nCol();
	const char * get_colName(int);
	const char * get_rowName(int);
	std::vector<const char *> get_colNames();
	std::vector<const char *> get_rowNames();
	void get_column(int, std::vector<double>&);			
};	

inline ExprMat::ExprMat(SEXP d):data(d){
	SEXP tmp;	
	dim=INTEGER(getAttrib(d, R_DimSymbol));
	
	PROTECT(tmp=VECTOR_ELT(getAttrib(data, R_DimNamesSymbol),0));
	for(int i=0;i<dim[0];++i){
		rownames.push_back( CHAR(STRING_ELT(tmp,i)) );
	}
	UNPROTECT(1);

	PROTECT(tmp=VECTOR_ELT(getAttrib(data, R_DimNamesSymbol),1));
	for(int i=0;i<length(tmp);++i){
		colnames.push_back( CHAR(STRING_ELT(tmp,i)) );
	}
	UNPROTECT(1);
}
 
inline  int ExprMat::get_nRow(){
	return(dim[0]);
}

inline int ExprMat::get_nCol(){
	return(dim[1]);
}

inline const char * ExprMat::get_colName(int index){//0 based index
	return( colnames[index] );
}

inline const char * ExprMat::get_rowName(int index){//0 based index
	return( rownames[index] );
}

inline std::vector <const char *> ExprMat::get_colNames(){//0 based index
	return( colnames );
}

inline std::vector <const char *> ExprMat::get_rowNames(){//0 based index
	return( rownames );
}

void ExprMat::get_column(int index, std::vector<double>& colmn){//0 based index
	double * dd =REAL(data);
	for(int i=(0+(dim[0]*index));i<(dim[0]*(index+1));++i){
		colmn.push_back(dd[i]);
	}
}


////////////////////////////////////////////////////////////////////////////////
 

class ProbeSets{
	private:	
	std::vector< const char * > prNames;
	std::vector< std::vector<int> > prIndex;
	std::vector< std::vector<int> > refIndex;
	public:	
	ProbeSets(SEXP,SEXP);
	int get_prSetCount();
	const char * get_prSetName(int);//args: index  0-based
	void get_prIndex(int, std::vector<int>&);//args: index  0-based	
	void get_refIndex(int, std::vector<int>&);//args: index  0-based		
};	

inline ProbeSets::ProbeSets(SEXP index, SEXP index_ref){
	prNames.reserve(length(index));
	prIndex.reserve(length(index));	
	refIndex.reserve(length(index_ref));	
	
	//prNames
	SEXP tmp;
	PROTECT(tmp=getAttrib(index, R_NamesSymbol)); // prNames
	prNames.reserve(length(tmp));
	for(int i=0;i<length(tmp);++i)
		prNames.push_back( CHAR(STRING_ELT(tmp,i)) );
	UNPROTECT(1);	

	int * tmpIndex;
	std::vector<int> tmpIndexVec;
	// prIndex
	for(int i=0;i<length(index);++i){ 
		tmpIndex=INTEGER(VECTOR_ELT(index,i));
		for(int j=0;j<length(VECTOR_ELT(index,i));++j){
			tmpIndexVec.push_back(tmpIndex[j]);
		}			
		prIndex.push_back(tmpIndexVec);
		std::vector<int>().swap(tmpIndexVec);		
	}
	// refIndex	
	for(int i=0;i<length(index_ref);++i){  	
		tmpIndex=INTEGER(VECTOR_ELT(index_ref,i));
		for(int j=0;j<length(VECTOR_ELT(index_ref,i));++j){
			tmpIndexVec.push_back(tmpIndex[j]);
		}			
		refIndex.push_back(tmpIndexVec);
		std::vector<int>().swap(tmpIndexVec);		
	}			
}

inline int ProbeSets::get_prSetCount(){ 
	return(prIndex.size()); 
}

inline const char * ProbeSets::get_prSetName(int index){ 
	return(prNames[index]); 
}

inline void ProbeSets::get_prIndex(int index, std::vector<int>& out){
	out=prIndex[index];
}	

inline void ProbeSets::get_refIndex(int index, std::vector<int>& out){
	out=refIndex[index];
}



