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
//      Pre declarations of matrix functions
//-----------------------------------------------------------------------
Matrix* new_matrix(unsigned int rows, unsigned int cols);
void free_matrix(Matrix* M); 
Matrix* multiply_matrix(Matrix* A, Matrix* B);
void print_matrix(Matrix* m);

//-----------------------------------------------------------------------
//      Pre declarations of test functions
//-----------------------------------------------------------------------
void test_multi();

//-----------------------------------------------------------------------
//      Main funciton
//-----------------------------------------------------------------------
int main(int argc, char *argv[]) {
    test_multi();
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
