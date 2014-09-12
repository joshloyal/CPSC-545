#include <stdio.h>
#include <stdlib.h>

struct Matrix {
    unsigned int rows;
    unsigned int cols;
    double* matrix; /* matrix stored in row major order */
};
#define M(m,x,y) m.matrix[ y + x*m.cols ]
typedef struct Matrix Matrix;

Matrix new_matrix(unsigned int rows, unsigned int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.matrix = (double*) malloc( rows*cols * sizeof(double) );

    return mat;
}

void print_matrix(Matrix m) {
    int i, j;
    for(i = 0; i < m.rows; i++ ) {
        printf("| ");
        for(j = 0; j < m.cols; j++) {
            printf("%f ", M(m, i, j));
        }
        printf("|\n");
    }
    printf("\n");
}

int main(int argc, char *argv[]) {
    int i,j;
    int size = 3;
    Matrix a;

    a = new_matrix(size,size);
    for(i = 0; i < size; i++ ) {
        for(j = 0; j < size; j++ ) {
            M(a,i,j) = (float)rand() / RAND_MAX;
        }
    }

    print_matrix(a);

    return 0;
}
