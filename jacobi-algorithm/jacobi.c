#include <stdio.h>
#include <stdlib.h>

/* Matrix Structure :
 *      Allows one to access elements via M->A[i][j] while hiding
 *      the implementation of storing the numbers in a single contigent
 *      row-major ordered array.
 */
struct Matrix {
    
    unsigned int rows, cols; // row and column dimensions of the matrix
    double* raw;             // the 'raw pointer' to the matrix elements
    double* *A;              // pointer to an array of pointers to rows
};
typedef struct Matrix Matrix;

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

void free_matrix(Matrix* M) {
    free(M->A);
    free(M->raw);
    free(M);
}

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

int main(int argc, char *argv[]) {
    int i,j;
    int size = 3;
    Matrix* a;

    a = new_matrix(size,size);
    for(i = 0; i < size; i++ ) {
        for(j = 0; j < size; j++ ) {
            a->A[i][j] = (float)rand() / RAND_MAX;
        }
    }

    print_matrix(a);
    
    free_matrix(a);

    return 0;
}
