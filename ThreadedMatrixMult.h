//
// Created by shepa on 3/10/2023.
//

#ifndef MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H
#define MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H

#include <pthread.h>
#include <thread>
#include <vector>
using namespace std;


struct SquareMatrix{
    int dim;
    int** data;    // points to a [dim x dim] square matrix
    SquareMatrix(int dim) {
        this->dim = dim;
        data = new int*[dim];
        for(int i = 0; i < dim; i++) {
            data[i] = new int[dim];
            for(int j = 0; j < dim; j++) {
                data[i][j] = 0;
            }
        }
    }
};

struct MatrixThreadParams {
    SquareMatrix* A;
    SquareMatrix* B;
    SquareMatrix* R;
    int rowOffset, columnOffset, dim;
};

void * BruteForceSquareMatrixMultiplication(void *package) {
    MatrixThreadParams* matrixPackage = (MatrixThreadParams*)package;

    //SquareMatrix* resM = new SquareMatrix{matrixPackage->A->dim, new int*[matrixPackage->A->dim]};

    for(int i = 0; i < matrixPackage->rowOffset; i++) {
        for(int j = 0; j < matrixPackage->columnOffset; j++) {
            matrixPackage->R->data[i][j] = 0;
            for(int k = 0; k < m->dim; k++) {
                resM->data[i][j] += matrixPackage->A->data[i][k] * matrixPackage->B->data[k][j];
            }
        }
    }
    return (void *)resM;
}
SquareMatrix* BruteForce(const SquareMatrix& A, const SquareMatrix& B) {
    SquareMatrix* m = new SquareMatrix(A.dim);
    int sum;
    for(int i = 0; i < m->dim; i++) {
        m->data[i] = new int[m->dim];
        for(int j = 0; j < m->dim; j++) {
            m->data[i][j] = 0;
            sum = 0;
            for(int k = 0; k < m->dim; k++) {
                sum += A.data[i][k] * B.data[k][j];
            }
            m->data[i][j] = sum;
        }
    }
    return m;
}
SquareMatrix* ThreadedDivideAndConquer(const SquareMatrix& A, const SquareMatrix& B) {

    SquareMatrix* resM = new SquareMatrix(A.dim);

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

        // spawn threads
        vector<thread> threads;

        // spawn 4 threads
        threads.emplace_back([&C11 = resM, &A11, &B11]() {C11 = ThreadedDivideAndConquer(A11, B11);});
        threads.emplace_back([&C12 = resM, &A12, &B12]() {C12 = ThreadedDivideAndConquer(A12, B12);});
        threads.emplace_back([&C21 = resM, &A21, &B21]() {C21 = ThreadedDivideAndConquer(A21, B21);});
        threads.emplace_back([&C22 = resM, &A22, &B22]() {C22 = ThreadedDivideAndConquer(A22, B22);});

        for(auto& t: threads) {
            t.join();
        }
    }
    return resM;

}
SquareMatrix* Strassen(const SquareMatrix& A, const SquareMatrix& B);
SquareMatrix* ThreadedStrassen(const SquareMatrix& A, const SquareMatrix& B) {

}



#endif//MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H
