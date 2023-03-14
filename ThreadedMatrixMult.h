//
// Created by shepa on 3/10/2023.
//

#ifndef MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H
#define MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H

#include <pthread.h>
#include <thread>
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
    const SquareMatrix* A;
    const SquareMatrix* B;
    SquareMatrix* R;
    int rowEnd, colEnd, dim;
    int rowStart, colStart;
    MatrixThreadParams(const SquareMatrix* A, const SquareMatrix* B, SquareMatrix* R, int rowEnd, int colEnd, int rowStart, int colStart, int dim) {
        this->A = A;
        this->B = B;
        this->R = R;
        this->rowEnd = rowEnd;
        this->colEnd = colEnd;
        this->rowStart = rowStart;
        this->colStart = colStart;
        this->dim = dim;
    }
};

SquareMatrix* matrixAddition(const SquareMatrix& A, const SquareMatrix& B) {
    SquareMatrix* C = new SquareMatrix(A.dim);

    for(int i = 0; i < A.dim; i++) {
        for(int j = 0; j < B.dim; j++) {
            C->data[i][j] = A.data[i][j] + B.data[i][j];
        }
    }
    return C;
}

SquareMatrix* matrixSubtraction(const SquareMatrix& A, const SquareMatrix& B) {
    SquareMatrix* C = new SquareMatrix(A.dim);

    for(int i = 0; i < A.dim; i++) {
        for(int j = 0; j < B.dim; j++) {
            C->data[i][j] = A.data[i][j] - B.data[i][j];
        }
    }
    return C;
}

void BruteForceSquareMatrixMultiplication(MatrixThreadParams* matrixParams) {
    //MatrixThreadParams* matrixParams = (MatrixThreadParams*)params;

    //SquareMatrix* resM = new SquareMatrix{matrixPackage->A->dim, new int*[matrixPackage->A->dim]};

    for(int i = matrixParams->rowStart; i < matrixParams->rowEnd; i++) {
        for(int j = matrixParams->colStart; j < matrixParams->colEnd; j++) {
            //matrixParams->R->data[i][j] = 0;
            int sum = 0;
            for(int k = 0; k < matrixParams->dim; k++) {
                sum += matrixParams->A->data[i][k] * matrixParams->B->data[k][j];
            }
            matrixParams->R->data[i][j] = sum;
        }
    }
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

        //DIVIDE AND CONQUER
        // spawn 4 threads
        //pthread_t id1, id2, id3, id4;
        int dim = resM->dim;
        MatrixThreadParams* R11Param = new MatrixThreadParams((SquareMatrix *) &A, (SquareMatrix*) &B, resM, mid, mid, 0, 0, dim);
        MatrixThreadParams* R12Param = new MatrixThreadParams((SquareMatrix *) &A, (SquareMatrix*) &B, resM, mid + mid, mid, mid, 0, dim);
        MatrixThreadParams* R21Param = new MatrixThreadParams((SquareMatrix *) &A, (SquareMatrix*) &B, resM, mid, mid + mid, 0, mid, dim);
        MatrixThreadParams* R22Param = new MatrixThreadParams((SquareMatrix *) &A, (SquareMatrix*) &B, resM, mid + mid, mid + mid, mid, mid, dim);

        /*
        cout << "Thread 1" << endl;
        pthread_create(&id1, NULL, &BruteForceSquareMatrixMultiplication, (void *)R11Param);

        cout << "Thread 2" << endl;
        pthread_create(&id2, NULL, &BruteForceSquareMatrixMultiplication, (void *)R12Param);

        cout << "Thread 3" << endl;
        pthread_create(&id3, NULL, &BruteForceSquareMatrixMultiplication, (void *)R21Param);

        cout << "Thread 4" << endl;
        pthread_create(&id4, NULL, &BruteForceSquareMatrixMultiplication, (void *)R22Param);
         */


        thread thread1(BruteForceSquareMatrixMultiplication, R11Param);
        thread thread2(BruteForceSquareMatrixMultiplication, R12Param);
        thread thread3(BruteForceSquareMatrixMultiplication, R21Param);
        thread thread4(BruteForceSquareMatrixMultiplication, R22Param);

        thread1.join();
        thread2.join();
        thread3.join();
        thread4.join();

    }
    return resM;

}
SquareMatrix* Strassen(const SquareMatrix& A, const SquareMatrix& B) {

    int mid = A.dim/2;
    SquareMatrix* A11 = new SquareMatrix(mid);
    SquareMatrix* A12 = new SquareMatrix(mid);
    SquareMatrix* A21 = new SquareMatrix(mid);
    SquareMatrix* A22 = new SquareMatrix(mid);
    SquareMatrix* B11 = new SquareMatrix(mid);
    SquareMatrix* B12 = new SquareMatrix(mid);
    SquareMatrix* B21 = new SquareMatrix(mid);
    SquareMatrix* B22 = new SquareMatrix(mid);

    // fill matrices
    for(int i = 0; i < mid; i++) {
        for(int j = 0; j < mid; j++) {
            A11->data[i][j] = A.data[i][j];
            A12->data[i][j] = A.data[i][j + mid];
            A21->data[i][j] = A.data[i + mid][j];
            A22->data[i][j] = A.data[i + mid][j + mid];

            // fill all B sub-matrices with data
            B11->data[i][j] = B.data[i][j];
            B12->data[i][j] = B.data[i][j + mid];
            B21->data[i][j] = B.data[i + mid][j];
            B22->data[i][j] = B.data[i + mid][j + mid];
        }
    }
    SquareMatrix* s1 = BruteForce(*A11, *matrixSubtraction(*B12, *B22));
    SquareMatrix* s2 = BruteForce(*matrixAddition(*A11, *A12), *B22);
    SquareMatrix* s3 = BruteForce(*matrixAddition(*A21, *A22), *B11);
    SquareMatrix* s4 = BruteForce(*A22, *matrixSubtraction(*B21, *B11));
    SquareMatrix* s5 = BruteForce(*matrixAddition(*A11, *A22), *matrixAddition(*B11, *B22));
    SquareMatrix* s6 = BruteForce(*matrixSubtraction(*A12, *A22), *matrixAddition(*B21, *B22));
    SquareMatrix* s7 = BruteForce(*matrixSubtraction(*A11, *A21), *matrixAddition(*B11, *B12));

    SquareMatrix* I = matrixAddition(*s5, *matrixAddition(*s6, *matrixSubtraction(*s4, *s2)));
    SquareMatrix* J = matrixAddition(*s1, *s2);
    SquareMatrix* K = matrixAddition(*s3, *s4);
    SquareMatrix* L = matrixSubtraction(*s1, *matrixSubtraction(*s7, *matrixAddition(*s3, *s5)));

    SquareMatrix* resM = new SquareMatrix(A.dim);

    for(int i = 0; i < mid; i++) {
        for(int j = 0; j < mid; j++) {
            resM->data[i][j] = I->data[i][j];
            resM->data[i][j + mid] = J->data[i][j];
            resM->data[i + mid][j] = K->data[i][j];
            resM->data[i + mid][j + mid] = L->data[i][j];
        }
    }
    return resM;

}
SquareMatrix* ThreadedStrassen(const SquareMatrix& A, const SquareMatrix& B) {

}




#endif//MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H
