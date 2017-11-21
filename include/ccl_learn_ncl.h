#ifndef __CCL_LEARN_NCL_H
#define __CCL_LEARN_NCL_H

#include <ccl_math.h>
#include <ccl_learn_alpha.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
    int dim_b;
    int dim_x;
    int dim_y;
    int dim_n;
    double* c;
    double s2;
    double* w;
}LEARN_NCL_MODEL;

typedef struct {
    gsl_matrix* HS;
    gsl_matrix* g;
    gsl_matrix* Y_T;
    gsl_matrix* Y_;
    gsl_matrix* H;
    gsl_matrix* BX_T;
    gsl_matrix* BX_;
    gsl_matrix* w_;
    gsl_matrix* pinvH1;
    gsl_vector* V;
    gsl_matrix* D;
    int    * idx;
}LEARN_MODEL_WS;


typedef struct {
    double* W;
    double* W_;
    gsl_matrix* J;
    gsl_vector* b_n;
    gsl_matrix* b_n_T;
    gsl_vector* u_n;
    gsl_matrix* u_n_T;
    gsl_matrix* BX_;
    gsl_matrix* Y_;
    gsl_vector* Wb;
    gsl_matrix* Wb_T;
    double c;
    double a;
    gsl_matrix* j_n;
    gsl_vector* j_n_flt;
    gsl_matrix* tmp2;
    double tmp;
} OBJ_WS;

typedef struct{
    int      dim_x;
    int      dim_n;
    int      dim_b;
    int      dim_y;
    int      r_ok;
    int      d_ok;
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
} SOLVE_NONLIN_WS;

int ccl_learn_ncl_model_alloc(LEARN_NCL_MODEL *model);
int ccl_learn_ncl_model_free(LEARN_NCL_MODEL *model);
void ccl_learn_ncl(const double * X, const double *Y, const int dim_x, const int dim_y, const int dim_n, const int dim_b, LEARN_NCL_MODEL *model);
void ccl_learn_model_dir(LEARN_NCL_MODEL *model, const double *BX, const double *Y);
void ccl_learn_model_ws_alloc(LEARN_NCL_MODEL *model,LEARN_MODEL_WS* ws);
void ccl_learn_model_ws_free(LEARN_MODEL_WS* ws);
void obj_ncl(const LEARN_NCL_MODEL *model, const double* W, const double*BX, const double*Y, double*fun,double* J);
void obj_ws_alloc(const LEARN_NCL_MODEL *model,OBJ_WS* ws);
void obj_ws_free(OBJ_WS* ws);
void ccl_lsqnonlin(const LEARN_NCL_MODEL* model,const  double* BX, const double*Y, const OPTION option,SOLVE_NONLIN_WS * lm_ws_param, double* W);
int ccl_solve_nonlin_ws_alloc(const LEARN_NCL_MODEL *model,SOLVE_NONLIN_WS * lm_ws);
int ccl_solve_nonlin_ws_free(SOLVE_NONLIN_WS * lm_ws);
void predict_ncl(const LEARN_NCL_MODEL* model, const double* BX, double* Unp);
int ccl_write_ncl_model(char* filename, LEARN_NCL_MODEL *model);
#ifdef __cplusplus
}
#endif
#endif

