#ifndef __CCL_MATH_H
#define __CCL_MATH_H

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort.h>

typedef struct {
    gsl_matrix *V;
    gsl_matrix *Sigma_pinv;
    gsl_matrix *U;
    gsl_matrix *A_pinv;
    gsl_matrix *A;
    gsl_matrix *_tmp_mat;
    gsl_vector *_tmp_vec;
    gsl_vector *u;
}MP_INV_WS;

void ccl_mat_add (double *A, const double *B, int row, int col);
void ccl_mat_sub(double *A, const double *B, int row, int col);
double ccl_vec_sum(double* vec, int size);
void mat_hotz_app ( double *A, int i,int j, const double *B,int k,int d,double* c);
void mat_vert_app (const double *A, int i, int j, const double *B, int k, int d, double * c);
void ccl_dot_product (const double *A, int i,int j, const double *B,int k,int d,double *C);
void ccl_mat_inv(const double *A, int i, int j, double *invA);
void ccl_MP_pinv(const double *A_in, int row, int col, double *invA);
void ccl_MP_inv_ws_alloc(MP_INV_WS *ws, int n, int m);
void ccl_MP_inv_ws_free(MP_INV_WS *ws);
void ccl_MP_pinv_test(const double *A_in, int row, int col, MP_INV_WS *ws,double *invA);

void linspace(double min, double max, double n, double *y);
void repmat(const double *mat, int mat_r, int col, int rows, int cols, double *A);
void repvec(const gsl_vector *vec, int rows, int cols, gsl_matrix * A);
void repvvec(const double *vec, int size, int rows, double *A);
void flt_mat(const double *mat, int row, int col, double *vec);
void ccl_mat_var(const double* data_in,int row,int col,int axis,double * var);
void print_mat(gsl_matrix *mat);
void print_mat_d(double *mat, int row, int col);
void print_mat_i(int * mat,int row,int col);
void vec_to_mat(const gsl_vector *vec, gsl_matrix *mat);
void nround(const double *n, int size, unsigned int c, double *ret);
void generate_kmeans_centres(const double * X, const int dim_x,const int dim_n, const int dim_b,double * centres);
void ccl_get_sub_mat_cols(const double * mat, const int row, const int col,const int * ind, int size, double * ret);
int ccl_mat_distance(const double *A,const int a_i,const int a_j,const double *B,const int b_i,const int b_j,double * D);
void ccl_mat_sum(const double *mat, const int i, const int j, const int axis, double * ret);
void ccl_mat_min(const double * mat,const int i,const int j,const int axis,double* val,int* indx);
void ccl_mat_mean(const double *mat,const int i, const int j,const int axis,double*val);
int ccl_find_index_int(const int *a, const int num_elements, const int operand, const int value,int* idx);
int ccl_find_index_double(const double *a, const int num_elements, const int operand, const double value,int* idx);
void ccl_mat_set_col(double * mat, int i, int j, int c, double* vec);
void ccl_gaussian_rbf(const double * X,int i, int j,const double *C,int k, int d,double s,double * BX);
void ccl_mat_transpose(const double* mat, int i, int j, double* mat_T);
int  ccl_any(double* vec, double epx, int size, int op);
void ccl_mat_reshape(const double *vec, int i, int j, double *mat);
#ifdef __cplusplus
}
#endif
#endif

