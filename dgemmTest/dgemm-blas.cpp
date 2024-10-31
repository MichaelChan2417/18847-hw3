#ifdef __APPLE__
#include <Accelerate/Accelerate.h> 
#else
#include <cblas.h>
#endif

const char* dgemm_desc = "BLAS dgemm.";
void square_dgemm( int n, double *A, double *B, double *C ) {
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A, n, B, n, 0.0, C, n);
}
