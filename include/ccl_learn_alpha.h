#ifndef __CCL_LEARN_ALPHA_H
#define __CCL_LEARN_ALPHA_H

#include <ccl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define NUM_CONSTRAINT 3
#ifdef __cplusplus
extern "C" {
#endif
typedef struct{
    int dim_b;
    int dim_r;
    int dim_x;
    int dim_u;
    int dim_t;
    int dim_k;
    int dim_n;
    double var;
    double nmse;
    double *w[NUM_CONSTRAINT];
    double * c;
    double  s2;
} LEARN_A_MODEL;

typedef struct{
    int MaxIter;
    double Tolfun;
    double Tolx;
    double Jacob;
} OPTION;

typedef struct{
    int      dim_x;
    int      dim_u;
    int      dim_n;
    int      dim_b;
    int      dim_k;
    double * xc;
    double * x;
    double * xf;
    double * epsx;
    double   epsf;
    double * r;
    double * J;
    double   S;
    double * A;
    double * v;
    double * D;
    double   Rlo;
    double   Rhi;
    double   l;
    double   lc;
    double * d;
    int      iter;
    int      r_ok;
    int      x_ok;
    double * xd;
    double * rd;
    double   Sd;
    double   dS;
    double   R;
    double   nu;
    double*  d_T;
    double*  J_T;
    double* tmp;
    double* rd_T;
    gsl_matrix* D_pinv;
    gsl_vector* A_d;
    gsl_matrix* A_inv;
    gsl_vector* A_inv_diag;
    gsl_vector* r_T;
} SOLVE_LM_WS;
int ccl_learn_alpha_model_alloc(LEARN_A_MODEL *model);
int ccl_learn_alpha_model_free(LEARN_A_MODEL *model);
int ccl_solve_lm_ws_alloc(const LEARN_A_MODEL *model,SOLVE_LM_WS * lm_ws);
int ccl_solve_lm_ws_free(SOLVE_LM_WS * lm_ws);
void ccl_learn_alpha(const double * Un,const double *X,const int dim_b,const int dim_r,const int dim_n,const int dim_x,const int dim_u,LEARN_A_MODEL optimal);
void search_learn_alpha(const double *BX,const double *RnUn, LEARN_A_MODEL* model);
void obj_AUn (const LEARN_A_MODEL* model, const double* W, const double* BX,const double * RnUn,double* fun_out);
void ccl_get_unit_vector_from_matrix(const double *theta, int dim_n, int dim_t, double *alpha);
void ccl_solve_lm(const LEARN_A_MODEL* model,const  double* RnUn,const  double* BX, const OPTION option,SOLVE_LM_WS * lm_ws_param, double* W);
void findjac(const LEARN_A_MODEL* model, const int dim_x,const double* BX, const double * RnUn,const double *y,const double*x,double epsx,double* J);
void ccl_get_rotation_matrix(const double*theta,const double*currentRn,const LEARN_A_MODEL* model,const int alpha_id,double*Rn);
void ccl_make_given_matrix(const double theta,int i,int j,int dim,double*G);
void predict_proj_alpha(double* x, LEARN_A_MODEL* model,double* centres,double variance,double* Iu, double*A);
#ifdef __cplusplus
}
#endif
#endif

