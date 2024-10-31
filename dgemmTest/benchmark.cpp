#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>

#include <cfloat>
#include <cmath>

#include <chrono>
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif
using namespace std;
/*Your function must have the following signature: */


extern void square_dgemm( int M, double *A, double *B, double *C );

/* Helper functions */

inline double read_timer( )
{
    static auto start = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
    return elapsed.count();
}

void fill( double *p, int n )
{
    for( int i = 0; i < n; i++ )
        p[i] = 2 * drand48( ) - 1;
}

void absolute_value( double *p, int n )
{
    for( int i = 0; i < n; i++ )
        p[i] = fabs( p[i] );
}


/* The benchmarking program */

int main( int argc, char **argv )
{

    /* These sizes should highlight performance dips at multiples of certain powers-of-two */
    unsigned int test_sizes[] = {
      31, 32, 33, 63, 64, 65, 95, 96, 97, 127, 128, 129, 191, 192, 229, 255, 256, 257,
      319, 320, 321, 417, 479, 480, 511, 512, 513
    };
    unsigned int count=sizeof(test_sizes)/sizeof(test_sizes[0]);
    unsigned int MaxN = test_sizes[count-1];

    double *A = (double*) malloc( MaxN * MaxN * sizeof(double) );
    double *B = (double*) malloc( MaxN * MaxN * sizeof(double) );
    double *C = (double*) malloc( MaxN * MaxN * sizeof(double) );
    fill( A, MaxN * MaxN );
    fill( B, MaxN * MaxN );
 
    /*For each test size*/
    for( unsigned int isize = 0; isize < count; isize++ )
      {
	/*Create and fill 3 random matrices A,B,C*/
        int n = test_sizes[isize];


        
        /*  measure Mflop/s rate; time a sufficiently long sequence of calls to eliminate noise*/
        double Mflop_s, seconds = -1.0;
        for( int n_iterations = 1; seconds < 0.1; n_iterations *= 2 ) 
        {
            /* warm-up */
            square_dgemm( n, A, B, C );
            
            /*  measure time */
            seconds = read_timer( );
            for( int i = 0; i < n_iterations; i++ )
                square_dgemm( n, A, B, C );
            seconds = read_timer( ) - seconds;
           
            /*  compute Mflop/s rate */
            Mflop_s = 2e-6 * n_iterations * n * n * n / seconds;
        }
        //printf ("%d %g \n", n, Mflop_s);
        cout << "n " << n << ", MFlop/sec = " << Mflop_s << endl;
        /*  Ensure that error does not exceed the theoretical error bound */
		
		/* Set initial C to 0 and do matrix multiply of A*B */
        memset( C, 0, sizeof( double ) * n * n );
        square_dgemm( n, A, B, C );
		/*Subtract A*B from C using standard dgemm (note that this should be 0 to within machine roundoff)*/
        cblas_dgemm( CblasColMajor,CblasNoTrans,CblasNoTrans, n,n,n, -1, A,n, B,n, 1, C,n );
		/*Subtract the maximum allowed roundoff from each element of C*/
        absolute_value( A, n * n );
        absolute_value( B, n * n );
        absolute_value( C, n * n );
        cblas_dgemm( CblasColMajor,CblasNoTrans,CblasNoTrans, n,n,n, -3.0*DBL_EPSILON*n, A,n, B,n, 1, C,n );
		/*After this test if any element in C is still positive something went wrong in square_dgemm*/
        for( int i = 0; i < n * n; i++ )
          if( C[i] > 0 )
            {
                printf( "FAILURE: error in matrix multiply exceeds an acceptable margin\n" );
                exit(-1);
            }

		/*Deallocate memory*/

    }
    free( C );
    free( B );
    free( A );
    return 0;
}
