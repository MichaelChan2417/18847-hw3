#ifdef __APPLE__
#include <Accelerate/Accelerate.h> 
#else
#include <cblas.h>
#endif

const char* dgemm_desc = "BLAS dgemm.";
//void square_dgemm( int n, double *A, double *B, double *C )
//{
//}
