#ifndef __CCL_LEARN_NHAT_H
#define __CCL_LEARN_NHAT_H

#include <ccl_math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define NUM_SEARCH 27000 // NUM_SEARCH = num_theta^{dim_t};
#define NUM_CONSTRAINT 4 //NUM_CONSTRAINT = dim_u
#ifdef __cplusplus
extern "C" {
#endif
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

typedef struct {
    int dim_u;
    int dim_n;
    int num_theta;
    int dim_t;
    int dim_s;
    double epsilon;
    double *min_theta;
    double *max_theta;
    double *list;
    double *I_u;
    double *theta[NUM_SEARCH];
    double *alpha[NUM_SEARCH];
} NHAT_search;

typedef struct{
    int dim_u;
    int dim_n;
    int dim_t;
    int dim_c;

    double *theta;
    double *alpha;
    double *P;
    double variance;
    double umse_j;
    double nmse_j;
}NHAT_Model;

typedef struct{
    NHAT_Model model[NUM_CONSTRAINT];
}NHAT_result;

int  init_search_param(NHAT_search *search, int dim_u, int dim_n, int num_theta);
int  nhat_mem_alloc_search(NHAT_search *search);
int  nhat_mem_free_search(NHAT_search *search);

int  nhat_mem_alloc_model(NHAT_Model *model,const NHAT_search *search);
int  nhat_duplicate_model(NHAT_Model *dest, const NHAT_Model * src);
int  nhat_mem_free_model(NHAT_Model *model);

void get_unit_vector(const double* theta, int dim_t,double *alpha);
void generate_search_space(NHAT_search *search);
void search_first_alpha(const double *Vn,  const double *Un, NHAT_Model *model, const NHAT_search *search, double *stats);
void search_alpha_nhat(const double *Vn,  const double *Un, NHAT_Model *model, const NHAT_search *search, NHAT_Model *model_out,double *stats);
void learn_nhat(const double *Un, const int dim_u, const int dim_n, NHAT_Model *optimal);
void calclate_N(double * N, const double *A, int row, int col);
int  nhat_mem_alloc_result(NHAT_Model *model, const NHAT_search search, int dim_c);
int  nhat_mem_free_result(NHAT_result *result, int n_models);
int  nhat_mem_alloc_optimal(NHAT_Model *optimal,const NHAT_search search,int dim_c);
int  nhat_mem_free_optimal(NHAT_Model *optimal);
#ifdef __cplusplus
}
#endif
#endif

