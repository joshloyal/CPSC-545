#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//-----------------------------------------------------------------------
//  Matrix Structure :
//      Allows one to access elements via M->A[i][j] while hiding
//      the implementation of storing the numbers in a single contigent
//      row-major ordered array.
//-----------------------------------------------------------------------
typedef struct matrix {    
    unsigned int rows, cols; // row and column dimensions of the matrix
    double* raw;             // the 'raw pointer' to the matrix elements
    double* *A;              // pointer to an array of pointers to rows
} matrix;

//-----------------------------------------------------------------------
//  Jacobi Structure :
//      Hold jacobi rotation parameters
//-----------------------------------------------------------------------
typedef struct jacobi_param {
    double t;   // tangent
    double c;   // cosine 
    double s;   // sine
} jacobi_param;

//-----------------------------------------------------------------------
//  Singular Value Structure :
//      Holds singular values and their index in U,V 
//      (note: we re-order the singular values by size, so
//             U,V need to be re-ordered accordingly)
//-----------------------------------------------------------------------
typedef struct singular_value {
    unsigned int index;
    double value;
} singular_value;

//-----------------------------------------------------------------------
//      Pre declarations of matrix functions
//-----------------------------------------------------------------------
matrix* new_matrix(unsigned int rows, unsigned int cols, double* a);
void free_matrix(matrix* M); 
void print_matrix(matrix* m);

//-----------------------------------------------------------------------
//      Pre declarations of utility functions
//-----------------------------------------------------------------------
double sign(double val);
matrix* eye(unsigned int n);
double sum_entry(matrix* M, unsigned int i, unsigned int j);

//-----------------------------------------------------------------------
//      Pre declarations of jacobi-algorithm functions
//-----------------------------------------------------------------------
void jacobi(double* a, int n, double* s, double* u, double* v);
jacobi_param jacobi_parameters(double a_pp, double a_pq, double a_qq);
void rotate(matrix* L, unsigned int p, unsigned int q, matrix* R);

//-----------------------------------------------------------------------
//      Pre declarations of test functions
//-----------------------------------------------------------------------
void test_jac_rot();
void test_jacobi();

//-----------------------------------------------------------------------
//      Main funciton
//-----------------------------------------------------------------------
int main(int argc, char *argv[]) {
    test_jac_rot();
    test_jacobi();
    return 0;
}

//-----------------------------------------------------------------------
//      Implementation of matrix funciton
//-----------------------------------------------------------------------
matrix* new_matrix(unsigned int rows, unsigned int cols, double* a) {
    matrix* M;
    int i;

    M = malloc(sizeof(matrix));
    if( M == NULL)
        return M;

    M->rows = rows;
    M->cols = cols;

    // make the raw array with malloc
    if( a == NULL ) {
        M->raw = malloc( rows*cols * sizeof(double) );
    }
    else {
        M->raw = a;
    } 
    if (M->raw == NULL) {
        free(M);
        return NULL;
    }

    // create and initialize the array of row pointers
    M->A = malloc( rows * sizeof(double*) );
    if( M->A == NULL) {
        free(M->raw);
        free(M);
        return NULL;
    }
    for(i = 0; i < rows; i++) {
        // get a pointer to the start of row i
        M->A[i] = &(M->raw[i*cols]);
    }

    return M;
}

//-----------------------------------------------------------------------
void free_matrix(matrix* M) {
    free(M->A);
    free(M->raw);
    free(M);
}

//-----------------------------------------------------------------------
jacobi_param jacobi_parameters(double a_pp, double a_pq, double a_qq) {
    jacobi_param jac = {0., 0., 0.};

    double t = (a_pp - a_qq) / (2. * a_pq);
    jac.t = sign(t) / ( fabs(t) + sqrt(1 + t*t) );
    jac.c = 1 / sqrt(1 + jac.t*jac.t);
    jac.s = (jac.c * jac.t);

    return jac;
}

//-----------------------------------------------------------------------
void rotate(matrix* L, unsigned int p, unsigned int q, matrix* R) {
    // compute the jacobi rotation parameters
    double u_pp = sum_entry(L, p, p);
    double u_pq = sum_entry(L, p, q);
    double u_qq = sum_entry(L, q, q);
    jacobi_param jac = jacobi_parameters(u_pp, u_pq, u_qq);

    // update columns p and q of M,U
    int k; double temp;
    for(k = 0; k < L->rows; k++)
    {
        temp = L->A[k][p];
        L->A[k][p] = jac.s * L->A[k][q] + jac.c * temp;
        L->A[k][q] = jac.c * L->A[k][q] - jac.s * temp;

        temp = R->A[k][p];
        R->A[k][p] = jac.s * R->A[k][q] + jac.c * temp;
        R->A[k][q] = jac.c * R->A[k][q] - jac.s * temp;

    }
}

//-----------------------------------------------------------------------
double sign(double val) {
    if ( val < 0 ) {
        return -1;
    }
    return 1;
}

//-----------------------------------------------------------------------
matrix* eye(unsigned int n) {
    matrix* I = new_matrix(n, n, NULL);
    int i, j;
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            if (i != j ) {
                I->A[i][j] = 0.;
            }
            else {
                I->A[i][j] = 1.;
            }
        }
    }

    return I;
}

//-----------------------------------------------------------------------
double sum_entry(matrix* M, unsigned int i, unsigned int j) {
    int k;
    double sum = 0;
    for(k = 0; k < M->rows; k++) 
    {
        sum += M->A[k][i] * M->A[k][j];
    }
    return sum;
}

//-----------------------------------------------------------------------
void print_matrix(matrix* m) {
    int i, j;
    for(i = 0; i < m->rows; i++ ) {
        printf("| ");
        for(j = 0; j < m->cols; j++) {
            printf("%f ", m->A[i][j]);
        }
        printf("|\n");
    }
    printf("\n");
}

//-----------------------------------------------------------------------
/**
    Implements the Jacobi algorithm for the SVD of a real matrix
    
    input parameters:
    @param a the nxn matrix to be diagonalized, given as an nxn array;
             will ge destroyed by the subroutine.
    @param n the dimensionality of the matrix

    output parameters:
    @param s the spectrum of the matrix a, in order of decreasing values,
             given as an array of length n
    @param u the left singular values of the matrix a, as an nxn array
    @param v the right singular values of the matrix a, given as an nxn array
*/
void jacobi(double* a, int n, double* s, double* u, double* v) {
    // paramaters
    int i, j;
    double epsilon = 1e-15, comparison;
    int iter_max = 100000, iters = 0, count = 1;

    matrix* L = new_matrix(n, n, a); // -> matrix columns are left singular vectors
    matrix* R = eye(n); // -> matrix whose columns are right singular vectors
    
    // begin performing jacobi rotations
    while( (count != 0) && (iters < iter_max) ) { 

        // loop through all pairs i < j
        count = 0;
        for (i = 0; i < (n - 1); i++) 
        {
            for( j = i + 1; j < n; j++)  
            {
                comparison = epsilon * sqrt(sum_entry(L,i,i) * sum_entry(L,j,j));
                if( fabs(sum_entry(L,i,j)) > comparison) 
                {
                    rotate(L, i, j, R);
                    count++;
                }
            }
        }
        iters++;
    }
    
    // print information about the algorithm
    if (iters == iter_max) 
    {
        printf("\nJacobi's Method failed to converge. Terminating\n");
        return;
    }
    else
    {
        printf("\nJacobi's method converged in %i iterations.\n", iters);
    }

    // extract the singular values:
    //      left svds: norms of columns of L
    //      right svds: norms of columns of R
}



//-----------------------------------------------------------------------
//     Implementation of test functions 
//-----------------------------------------------------------------------
void test_jac_rot() {
    printf("//------------------------------------------------------------------\n");
    printf("Testing Jacobi rotation\n");
    printf("//------------------------------------------------------------------\n");
    jacobi_param jac = jacobi_parameters(2,4,3);
    printf("t = %f, s = %f, c = %f\n", jac.t, jac.s, jac.c);
}

void test_jacobi() {
    printf("//------------------------------------------------------------------\n");
    printf("Testing Jacobi Algorithm\n");
    printf("//------------------------------------------------------------------\n");
    unsigned int n = 3.;
    double* m = malloc( sizeof(double)*n*n );
    m[0] = 1.;
    m[1] = 3.;
    m[2] = 2.;
    m[3] = 5.;
    m[4] = 6.;
    m[5] = 4.;
    m[6] = 7.;
    m[7] = 8.;
    m[8] = 9.;
    matrix* M = new_matrix(n,n,m);

    printf("matrix M:\n");
    print_matrix(M);

    double* s = malloc( sizeof(double)*3 );
    double* u = malloc( sizeof(double)*9 );
    double* v = malloc( sizeof(double)*9 );
    jacobi(m, n, s, u, v);
}
