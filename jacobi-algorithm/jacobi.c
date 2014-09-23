#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//-----------------------------------------------------------------------
//  Matrix Structure :
//      Allows one to access elements via M->A[i][j] while hiding
//      the implementation of storing the numbers in a single contigent
//      row-major ordered array.
//-----------------------------------------------------------------------
struct Matrix {
    
    unsigned int rows, cols; // row and column dimensions of the matrix
    double* raw;             // the 'raw pointer' to the matrix elements
    double* *A;              // pointer to an array of pointers to rows
};
typedef struct Matrix Matrix;

//-----------------------------------------------------------------------
//  Jacobi Structure :
//      Hold jacobi rotation parameters
//-----------------------------------------------------------------------
struct Jacobi {
    double t;   // tangent
    double c;   // cosine 
    double s;   // sine
};
typedef struct Jacobi Jacobi;


//-----------------------------------------------------------------------
//      Pre declarations of matrix functions
//-----------------------------------------------------------------------
Matrix* new_matrix(unsigned int rows, unsigned int cols);
void free_matrix(Matrix* M); 
Matrix* multiply_matrix(Matrix* A, Matrix* B);
Matrix* transpose(Matrix* M);
Jacobi jacobi_rotation(double a_pp, double a_pq, double a_qq);
void print_matrix(Matrix* m);

//-----------------------------------------------------------------------
//      Pre declarations of utility functions
//-----------------------------------------------------------------------
double sign(double val);

//-----------------------------------------------------------------------
//      Pre declarations of test functions
//-----------------------------------------------------------------------
void test_multi();
void test_trans();
void test_jac_rot();

//-----------------------------------------------------------------------
//      Main funciton
//-----------------------------------------------------------------------
int main(int argc, char *argv[]) {
    test_multi();
    test_trans();
    test_jac_rot();
    return 0;
}

//-----------------------------------------------------------------------
//      Implementation of matrix funciton
//-----------------------------------------------------------------------
Matrix* new_matrix(unsigned int rows, unsigned int cols) {
    Matrix* M;
    int i;

    M = malloc(sizeof(Matrix));
    if( M == NULL)
        return M;

    M->rows = rows;
    M->cols = cols;

    // make the raw array with malloc
    M->raw = malloc( rows*cols * sizeof(double) );
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
void free_matrix(Matrix* M) {
    free(M->A);
    free(M->raw);
    free(M);
}

//-----------------------------------------------------------------------
Matrix* multiply_matrix(Matrix* A, Matrix* B) {
    unsigned int n, m, p;
    double sum;
    Matrix* C;
    
    if( A->cols != B->rows ) {
        printf("Dimension mismatch when multiplying matrices!\n", stderr);
        return NULL;
    }

    n = A->rows;
    m = A->cols;
    p = B->cols; 
    C = new_matrix(n,p);
    if( C != NULL ) {
        unsigned int i,j,k; 
        for(i = 0; i < n; i++) {
            for(j = 0; j < p; j++) {
                sum = 0;
                for(k = 0; k < m; k++) {
                    sum = sum + A->A[i][k] * B->A[k][j];
                }
                C->A[i][j] = sum;
            }
        }
    }
    
    return C;
}

/* just invert coordinates... */
Matrix* transpose(Matrix* M) {
    unsigned int rows, cols;
    Matrix* Mt;

    rows = M->rows;
    cols = M->cols;
    Mt = new_matrix(cols, rows); // flip rows and cols

    if ( Mt != NULL ) {
        unsigned int i, j;
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                Mt->A[j][i] = M->A[i][j];
            }
        }
    }

    return Mt;
}

//-----------------------------------------------------------------------
Jacobi jacobi_rotation(double a_pp, double a_pq, double a_qq) {
    Jacobi jac = {0., 0., 0.};

    double t = (a_qq - a_pp) / (2. * a_pq);
    jac.t = sign(t) / ( fabs(t) + sqrt(1 + t*t) );
    jac.c = 1 / sqrt(1 + jac.t*jac.t);
    jac.s = jac.c * jac.t;

    return jac;
}

//-----------------------------------------------------------------------
double sign(double val) {
    if ( val < 0 ) {
        return -1;
    }
    return 1;
}

//-----------------------------------------------------------------------
void rotate(Matrix* M) {
}

//-----------------------------------------------------------------------
void print_matrix(Matrix* m) {
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
//     Implementation of test functions 
//-----------------------------------------------------------------------
void test_multi() {
    printf("//------------------------------------------------------------------\n");
    printf("Testing Matrix Multiplication\n");
    printf("//------------------------------------------------------------------\n");
    Matrix* A;
    Matrix* B;
    Matrix* C;

    A = new_matrix(2,2);
    B = new_matrix(2,3);

    A->A[0][0] = 1;
    A->A[0][1] = 2;
    A->A[1][0] = 3;
    A->A[1][1] = 4;

    B->A[0][0] = 5;
    B->A[0][1] = 6;
    B->A[0][2] = 7;
    B->A[1][0] = 8;
    B->A[1][1] = 9;
    B->A[1][2] = 10;

    C = multiply_matrix(A,B);

    print_matrix(C); 

    free_matrix(A); free_matrix(B); free_matrix(C);
}

void test_trans() {
    printf("//------------------------------------------------------------------\n");
    printf("Testing Transpose\n");
    printf("//------------------------------------------------------------------\n");
    
    Matrix* M;
    Matrix* Mt;

    M = new_matrix(3,3);
    
    M->A[0][0] = 5;
    M->A[0][1] = 6;
    M->A[0][2] = 7;
    M->A[1][0] = 8;
    M->A[1][1] = 9;
    M->A[1][2] = 10;
    M->A[2][0] = 10;
    M->A[2][1] = 11;
    M->A[2][2] = 12;

    Mt = transpose(M);
    print_matrix(M);
    print_matrix(Mt);

    free_matrix(M); free_matrix(Mt);
}

void test_jac_rot() {
    printf("//------------------------------------------------------------------\n");
    printf("Testing Jacobi rotation\n");
    printf("//------------------------------------------------------------------\n");
    Jacobi jac = jacobi_rotation(2,4,3);
    printf("t = %f, s = %f, c = %f\n", jac.t, jac.s, jac.c);
}
