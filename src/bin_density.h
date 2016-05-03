
#pragma once

#include <R.h>
#include <Rdefines.h>


#define max(a,b)  ((a < b) ?  (b) : (a))
#define printf Rprintf

int vect_max(int * cpos,int wwidth,int * orivec);
int median(int * cpos,int wwidth,int * orivec );
int mean(int * cpos,int wwidth,int * orivec);
void shrink(int * orivec, int * newvec, int orivecl, int window_count, int (*summarizep)(int *,int,int *));
void expand(int * orivec, int * newvec, int orivecl, int window_count);


double vect_max_dble(int * cpos,int wwidth,double * orivec);
double median_dble(int * cpos,int wwidth,double * orivec );
double mean_dble(int * cpos,int wwidth,double * orivec);
void shrink_dble(double * orivec, double * newvec, int orivecl, int window_count, double (*summarizep)(int *,int,double *));
void expand_dble(double * orivec, double * newvec, int orivecl, int window_count);

SEXP approx_window(SEXP window_count, SEXP score_list, SEXP method);
