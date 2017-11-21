#ifndef __CCL_LEARN_POLICY_H
#define __CCL_LEARN_POLICY_H

#include <ccl_math.h>
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
#define NUM_CENTRES 20
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
}LEARN_MODEL_PI_WS;

typedef struct{
    int      dim_y;
    int      dim_x;
    int      dim_n;
    int      dim_b;
    double * w;
}LEARN_MODEL_PI;

typedef struct{
    int      dim_y;
    int      dim_x;
    int      dim_n;
    int      dim_b;
    int      dim_phi;
    double * c;
    double   s2;
    double * w[NUM_CENTRES];
}LEARN_MODEL_LW_PI;

typedef struct {
    gsl_vector* g;
    gsl_matrix* Y_N;
    gsl_vector* YN_vec;
    gsl_matrix* Y_Phit;
    gsl_matrix* ones;
    gsl_matrix* Y_;
    gsl_matrix* H;
    gsl_matrix* Phi;
    gsl_vector* Phi_vec;
    gsl_matrix* Phi_vec_T;
    gsl_matrix* YN_Phit;
    gsl_vector* YN_Phi_vec;
    gsl_matrix* YN_Phi_vec_T;
    gsl_matrix* vv;
    gsl_matrix* WX_;
    gsl_vector* WX_row;
    gsl_matrix* WPhi;
    gsl_matrix* WPhi_T;
    gsl_matrix* pinvH1;
    gsl_vector* V;
    gsl_vector* r;
    gsl_matrix* r_rep;
    gsl_matrix* D;
    int    * idx;
    gsl_vector* w_vec;
    gsl_matrix* w_;
    gsl_matrix* w_T;
    gsl_matrix* w[NUM_CENTRES];
}LEARN_MODEL_LW_PI_WS;
void ccl_learn_policy_pi(LEARN_MODEL_PI *model, const double *BX, const double *Y);
void ccl_learn_model_pi_ws_alloc(LEARN_MODEL_PI *model,LEARN_MODEL_PI_WS* ws);
void ccl_learn_model_pi_ws_free(LEARN_MODEL_PI_WS* ws);
void ccl_learn_policy_lw_pi(LEARN_MODEL_LW_PI *model, const double *WX, const double *X, const double *Y);
void ccl_learn_policy_lw_pi_model_alloc(LEARN_MODEL_LW_PI *model);
void ccl_learn_policy_lw_pi_model_free(LEARN_MODEL_LW_PI *model);
void ccl_learn_model_lw_pi_ws_alloc(LEARN_MODEL_LW_PI *model,LEARN_MODEL_LW_PI_WS* ws);
void ccl_learn_model_lw_pi_ws_free(LEARN_MODEL_LW_PI_WS* ws);
void predict_linear(const double* X, const double* centres,const double variance,const LEARN_MODEL_PI *model,double* Yp);
void predict_local_linear(const double* X, const double* centres,const double variance,const LEARN_MODEL_LW_PI *model,double* Yp);
int ccl_read_data_from_file(char* filename, int dim_x, int dim_n, double* mat);
int ccl_write_lwmodel_to_file(char* filename, LEARN_MODEL_LW_PI* model);
#ifdef __cplusplus
}
#endif
#endif

