#ifndef __GET_NMSE_H
#define __GET_NMSE_H

#include <ccl_config.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
     CCL_GAUSSIAN_KERNAL, CCL_BISQUARE_KERNAL
} CCL_Kernel;

void ccl_get_nmse(const int *a,const int * b);

#ifdef __cplusplus
}
#endif
#endif

