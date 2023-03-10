//
// Created by shepa on 3/10/2023.
//

#ifndef MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H
#define MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H

#include <thread>
using namespace std;


struct SquareMatrix{
    int dim;

    int** data;    // points to a [dim x dim] square matrix
};

struct Package {
    SquareMatrix matrix1;
    SquareMatrix matrix2;
};

void * BruteForceSquareMatrixMultiplication(void *package) {
    Package* matrixPackage = (Package*)package;

    SquareMatrix* resM = new SquareMatrix{matrixPackage->matrix1.dim, new int*[matrixPackage->matrix1.dim]};

    for(int i = 0; i < resM->dim; i++) {
        resM->data[i] = new int[resM->dim];
        for(int j = 0; j < resM->dim; j++) {
            resM->data[i][j] = 0;
            for(int k = 0; k < resM->dim; k++) {
                resM->data[i][j] += matrixPackage->matrix1.data[i][k] * matrixPackage->matrix2.data[k][j];
            }
        }
    }
    return (void *)resM;
}
SquareMatrix* BruteForce(const SquareMatrix& A, const SquareMatrix& B) {
    SquareMatrix* m = new SquareMatrix;
    m->dim = A.dim;
    m->data = new int*[m->dim];

    for(int i = 0; i < m->dim; i++) {
        m->data[i] = new int[m->dim];
        for(int j = 0; j < m->dim; j++) {
            m->data[i][j] = 0;
            for(int k = 0; k < m->dim; k++) {
                m->data[i][j] += A.data[i][k] * B.data[k][j];
            }
        }
    }
    return m;
}
SquareMatrix* ThreadedDivideAndConquer(const SquareMatrix& A, const SquareMatrix& B) {

    SquareMatrix* resM = new SquareMatrix;
    resM->dim = A.dim;
    if(A.dim == 1 && B.dim == 1) {
        resM->data[0][0] = A.data[0][0] * B.data[0][0];
    } else {
        int mid = A.dim / 2;

        SquareMatrix A11 = {mid, new int*[mid]};
        SquareMatrix A12 = {mid, new int*[mid]};
        SquareMatrix A21 = {mid, new int*[mid]};
        SquareMatrix A22 = {mid, new int*[mid]};
        SquareMatrix B11 = {mid, new int*[mid]};
        SquareMatrix B12 = {mid, new int*[mid]};
        SquareMatrix B21 = {mid, new int*[mid]};
        SquareMatrix B22 = {mid, new int*[mid]};

        // fill in matrices
        for(int i = 0; i < mid; i++) {
            A11.data[i] = new int[mid];
            A12.data[i] = new int[mid];
            A21.data[i] = new int[mid];
            A22.data[i] = new int[mid];
            B11.data[i] = new int[mid];
            B12.data[i] = new int[mid];
            B21.data[i] = new int[mid];
            B22.data[i] = new int[mid];

            for(int j = 0; j < mid; j++) {
                // fill all A sub-matrices with data
                A11.data[i][j] = A.data[i][j];
                A12.data[i][j] = A.data[i][j + mid];
                A21.data[i][j] = A.data[i + mid][j];
                A22.data[i][j] = A.data[i + mid][j + mid];

                // fill all B sub-matrices with data
                B11.data[i][j] = B.data[i][j];
                B12.data[i][j] = B.data[i][j + mid];
                B21.data[i][j] = B.data[i + mid][j];
                B22.data[i][j] = B.data[i + mid][j + mid];
            }
        }

        // spawn 4 threads
        SquareMatrix* res;
        thread thread1(res = (SquareMatrix*) BruteForceSquareMatrixMultiplication, new Package{A11, B11});
        thread1.join();

    }
    return resM;

}
SquareMatrix* Strassen(const SquareMatrix& A, const SquareMatrix& B);
SquareMatrix* ThreadedStrassen(const SquareMatrix& A, const SquareMatrix& B);



#endif//MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H
