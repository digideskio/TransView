
/*
 * bin_density.c
 *
 *  Created on: February 8, 2013
 *      Author: Julius Muller
 */

/**
 *  @file	bin_density.c
 *  @brief	Convenience functions to provide fast binning of the density vectors
 */


#include <math.h>
#include <signal.h>
#include "bin_density.h"

/**
* @brief Computes the max
*
* @param cpos Pointer to start index in orivec
* @param wwidth width of window
* @param orivec Original vector
* @return Max value
* @details
* @note Nothing
* @todo Nothing
*/
int vect_max(int * cpos,int wwidth,int * orivec){
	int mpos;
	int result=orivec[*cpos];
	for(mpos=*cpos;*cpos<mpos+wwidth;(*cpos)++){
		result=max(result,orivec[*cpos]);
	}
	return(result);
}

/**
* @brief Computes the median
*
* @param cpos Pointer to start index in orivec
* @param wwidth width of window
* @param orivec Original vector
* @return Median value
* @details Adopted from http://en.wikiversity.org/wiki/C_source_code_to_find_the_median_and_mean
* @note Nothing
* @todo Nothing
*/
int median(int * cpos,int wwidth,int * orivec ) {
	int temp;
    int j,mpos;
    for(mpos=*cpos; *cpos<mpos+wwidth; (*cpos)++) {
        for(j=*cpos+1; j<mpos+wwidth; j++) {
            if(orivec[j] < orivec[*cpos]) {
                temp = orivec[*cpos];
                orivec[*cpos] = orivec[j];
                orivec[j] = temp;
            }
        }
    }

    if(!(wwidth%2)){
    	return((orivec[*cpos-(wwidth/2)-1] + orivec[*cpos-((wwidth-1)/2)-1]) / 2.0);
    } else {
        return orivec[*cpos-(wwidth/2)-1];
    }
}

/**
* @brief Computes the mean
*
* @param cpos Pointer to start index in orivec
* @param wwidth width of window
* @param orivec Pointer to original vector
* @return Median value
* @details Simple mean algorithm
* @note Nothing
* @todo Nothing
*/
int mean(int * cpos,int wwidth,int * orivec){
	int mpos,msum=0;
	for(mpos=*cpos;*cpos<mpos+wwidth;(*cpos)++)msum+=orivec[*cpos];
	return(msum/wwidth);
}


/**
* @brief Summarizes a longer vector into bins of approximately equal width
*
* @param orivec Pointer to original vector
* @param wwidth Pointer to empty vector of size wsize
* @param orivecl length of the original vector
* @param wsize Window width of the vector that will be returned
* @param summarizep Pointer to the appropriate function that summarizes the windows such as mean median...
* @return Modifies newvec
* @details Walks through the vector window by window and calls summarizep each time
* @note Nothing
* @todo Nothing
*/
void shrink(int * orivec, int * newvec, int orivecl, int window_count, int (*summarizep)(int *,int,int *)){
	int cpos=0,wwidth=ceil((double)orivecl/(double)window_count),nvc=0;
	for(;window_count>0;window_count--){
		if(orivecl<wwidth)wwidth=orivecl;
		orivecl-=wwidth;
		newvec[nvc++]=(*summarizep)(&cpos,wwidth,orivec);
		if(orivecl%window_count)wwidth=orivecl/(window_count-1);
	}
}

/**
* @brief Expands a vector into a longer vector of approximately equally sized bins.
*
* @param orivec Pointer to original vector
* @param wwidth Pointer to empty vector of size wsize
* @param orivecl length of the original vector
* @param wsize Window width of the vector that will be returned
* @return Modifies newvec
* @details Walks through the vector and expands window by window
* @note Nothing
* @todo Nothing
*/
void expand(int * orivec, int * newvec, int orivecl, int window_count){
	int wwidth=floor((double)window_count/(double)orivecl),cpos=0,wc=0,nvc=0;
	int rem=window_count%orivecl;
	for(;orivecl>0;orivecl--){
		for(wc=0;wc<wwidth;wc++){
			newvec[nvc++]=orivec[cpos];
			if(rem){
				newvec[nvc++]=orivec[cpos];
				--rem;
			}
		}
		++cpos;
	}
}



/**
* @brief Computes the max
*
* @param cpos Pointer to start index in orivec
* @param wwidth width of window
* @param orivec Original vector
* @return Max value
* @details
* @note Nothing
* @todo Nothing
*/
double vect_max_dble(int * cpos,int wwidth,double * orivec){
	int mpos;
	double result=orivec[*cpos];
	for(mpos=*cpos;*cpos<mpos+wwidth;(*cpos)++){
		result=max(result,orivec[*cpos]);
	}
	return(result);
}

/**
* @brief Computes the median
*
* @param cpos Pointer to start index in orivec
* @param wwidth width of window
* @param orivec Original vector
* @return Median value
* @details Adopted from http://en.wikiversity.org/wiki/C_source_code_to_find_the_median_and_mean
* @note Nothing
* @todo Nothing
*/
double median_dble(int * cpos,int wwidth,double * orivec ) {
    double temp;
    int j,mpos;
    for(mpos=*cpos; *cpos<mpos+wwidth; (*cpos)++) {
        for(j=*cpos+1; j<mpos+wwidth; j++) {
            if(orivec[j] < orivec[*cpos]) {
                temp = orivec[*cpos];
                orivec[*cpos] = orivec[j];
                orivec[j] = temp;
            }
        }
    }

    if(!(wwidth%2)){
        return((orivec[*cpos-(wwidth/2)-1] + orivec[*cpos-((wwidth-1)/2)-1]) / 2.0);
    } else {
        return orivec[*cpos-(wwidth/2)-1];
    }
}

/**
* @brief Computes the mean
*
* @param cpos Pointer to start index in orivec
* @param wwidth width of window
* @param orivec Pointer to original vector
* @return Median value
* @details Simple mean algorithm
* @note Nothing
* @todo Nothing
*/
double mean_dble(int * cpos,int wwidth,double * orivec){
	int mpos;
	double msum=0;
	for(mpos=*cpos;*cpos<mpos+wwidth;(*cpos)++)msum+=orivec[*cpos];
	return(msum/wwidth);
}


/**
* @brief Summarizes a longer vector into bins of approximately equal width
*
* @param orivec Pointer to original vector
* @param wwidth Pointer to empty vector of size wsize
* @param orivecl length of the original vector
* @param wsize Window width of the vector that will be returned
* @param summarizep Pointer to the appropriate function that summarizes the windows such as mean median...
* @return Modifies newvec
* @details Walks through the vector window by window and calls summarizep each time
* @note Nothing
* @todo Nothing
*/
void shrink_dble(double * orivec, double * newvec, int orivecl, int window_count, double (*summarizep)(int *,int,double *)){
	int cpos=0,wwidth=ceil((double)orivecl/(double)window_count),nvc=0;
	for(;window_count>0;window_count--){
		if(orivecl<wwidth)wwidth=orivecl;
		orivecl-=wwidth;
		newvec[nvc++]=(*summarizep)(&cpos,wwidth,orivec);
		if(orivecl%window_count)wwidth=orivecl/(window_count-1);
	}
}

/**
* @brief Expands a vector into a longer vector of approximately equally sized bins.
*
* @param orivec Pointer to original vector
* @param wwidth Pointer to empty vector of size wsize
* @param orivecl length of the original vector
* @param wsize Window width of the vector that will be returned
* @return Modifies newvec
* @details Walks through the vector and expands window by window
* @note Nothing
* @todo Nothing
*/
void expand_dble(double * orivec, double * newvec, int orivecl, int window_count){
	int wwidth=floor((double)window_count/(double)orivecl),cpos=0,wc=0,nvc=0;
	int rem=window_count%orivecl;
	for(;orivecl>0;orivecl--){
		for(wc=0;wc<wwidth;wc++){
			newvec[nvc++]=orivec[cpos];
			if(rem){
				newvec[nvc++]=orivec[cpos];
				--rem;
			}
		}
		++cpos;
	}
}



/**
* @brief Summarizes a list of vectors into a list of binned vectors of equal length. Each vector bin summarizes an approximately equal amount of values.
*
* @param method Charater array defining the method to be used for binning. Can be 'mean' 'median' or 'max'
* @param score_list List with numeric vectors
* @param window_size Window width of the vectors that will be returned
* @return List with updated vectors
* @details Walks through the vectors and calls shrink or expand to set vectors to equal widths
* @note Nothing
* @todo Nothing
*/
SEXP approx_window(SEXP window_count, SEXP score_list, SEXP method) {
	const char *methodn = STRING_VALUE(method);
	const int wsize=INTEGER_VALUE(window_count);

	SEXP lnames = getAttrib(score_list, R_NamesSymbol);
	SEXP ori_vec,new_vec,out_names,out_list;
	int elcount=0,elements=LENGTH(lnames),upc=0,olen;
	signal(SIGINT,SIG_DFL);
	PROTECT(lnames = AS_CHARACTER(lnames));upc++;
	PROTECT(out_list = allocVector(VECSXP, elements));upc++;
	PROTECT(out_names = allocVector(STRSXP,elements));upc++;

	//Select proper call back
	double (*summarizep)(int *,int,double *);
	if(!strcmp(methodn,"mean")){
		summarizep=mean_dble;
	}else if(!strcmp(methodn,"median")){
		summarizep=median_dble;
	}else if(!strcmp(methodn,"max")){
		summarizep=vect_max_dble;
	}else{
		error("%s not known",methodn);
		goto FINALIZE;
	}


	for(;elcount<elements;++elcount){
	  PROTECT(ori_vec=AS_NUMERIC(VECTOR_ELT(score_list, elcount)));
	  PROTECT(new_vec = NEW_NUMERIC(wsize));
	  olen=LENGTH(ori_vec);
	  double *ori_vecp= NUMERIC_POINTER(ori_vec);
	  double *new_vecp= NUMERIC_POINTER(new_vec);
	  SET_STRING_ELT(out_names,elcount,mkChar(CHAR(STRING_ELT(lnames, elcount))));
	  if(olen>wsize){
		  shrink_dble(ori_vecp,new_vecp,olen,wsize,summarizep);
		  SET_VECTOR_ELT(out_list, elcount, new_vec);
	  }else if(olen<wsize){
		  expand_dble(ori_vecp,new_vecp,olen,wsize);
		  SET_VECTOR_ELT(out_list, elcount, new_vec);
	  }else{
		  SET_VECTOR_ELT(out_list, elcount, ori_vec);
	  }

	  UNPROTECT(2);
	}
	setAttrib(out_list, R_NamesSymbol, out_names);

	FINALIZE:
	UNPROTECT(upc);
	return(out_list);
}


