#ifndef __CCL_LEARN_LAMBDA_H
#define __CCL_LEARN_LAMBDA_H

#include <../include/ccl_math.h>
#include <../include/ccl_learn_alpha.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define NUM_CONSTRAINT 2
#ifdef __cplusplus
extern "C" {
#endif
typedef void (*JACOBIAN)(const double*,const int,double*);
void Jacobian(const double* X,const int size,double* out);
int ccl_learn_lambda_model_alloc(LEARN_A_MODEL *model);
int ccl_learn_lambda_model_free(LEARN_A_MODEL *model);
void ccl_learn_lambda(const double * Un,const double *X,void (*J_func)(const double*,const int,double*),const int dim_b,const int dim_r,const int dim_n,const int dim_x,const int dim_u,LEARN_A_MODEL optimal);
void predict_proj_lambda(double* x, LEARN_A_MODEL model,void (*J_func)(const double*,const int,double*),double* centres,double variance,double* Iu, double*A);
int ccl_write_learn_lambda_model(char* filename, LEARN_A_MODEL *model);
#ifdef __cplusplus
}
#endif
#endif

