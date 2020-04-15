static void HandError( cudaError_t err, const char *file, int line   ) {
        if (err != cudaSuccess) {
                    printf( "\n%s in %s at line %d\n", \
                                             cudaGetErrorString( err   ), file, line );
                             exit( EXIT_FAILURE   );
                                 }
         }
 #define HRR( err   ) (HandError( err, __FILE__, __LINE__   ))
