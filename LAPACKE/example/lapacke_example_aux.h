#ifndef _LAPACKE_EXAMPLE_AUX_
#define _LAPACKE_EXAMPLE_AUX_


void print_matrix_rowmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm );
void print_matrix_colmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm );
void print_vector( char* desc, lapack_int n, lapack_int* vec );

void print_matrix_rowmajor_64( char* desc, int64_t m, int64_t n, double* mat, int64_t ldm );
void print_matrix_colmajor_64( char* desc, int64_t m, int64_t n, double* mat, int64_t ldm );
void print_vector_64( char* desc, int64_t n, int64_t* vec );

#endif /* _LAPACKE_EXAMPLE_AUX_*/
