/*
* Author: Shepard Berry
* Assignment Title: MatrixMultiplication
* Assignment Description: This program demonstrates 4 different matrix multiplication
 *                          algorithms.
* Due Date: 3/26/2023
* Date Created: 3/3/2023
* Date Last Modified: 3/26/2023
 */

#ifndef MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H
#define MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H

#include <pthread.h>
using namespace std;

struct SquareMatrix{

    int dim;
    int** data;    // points to a [dim x dim] square matrix

    /*
     * description: Default SquareMatrix constructor
     * return: none
     * precondition: none
     * postcondition: A SquareMatrix is constructed
    */
    SquareMatrix() {}

    /*
     * description: SquareMatrix constructor with dimension parameter.
     * return: none
     * precondition: Dimensions necessary for SquareMatrix are declared
     * postcondition: A SquareMatrix is constructed
    */
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

struct MatrixSection {
    // refers to a complete matrix, and knows its bounds
    SquareMatrix* M;
    int rowStart, colStart, dim;

    /*
     * description: Constructs a matrix section with an underlying matrix
     *                  and bounded dimensions
     * return: none
     * precondition: A matrix in memory is declared
     * postcondition: A MatrixSection is constructed
    */
    MatrixSection(SquareMatrix* M, int rowStart, int colStart, int dim) {
        this->M = M;
        this->rowStart = rowStart;
        this->colStart = colStart;
        this->dim = dim;
    }
};
struct MatrixThreadParams {
    const MatrixSection* A;
    const MatrixSection* B;
    MatrixSection* R;

    /*
     * description: Constructs a MatrixThreadParams object
     * return: none
     * precondition: 3 matrices are declared
     * postcondition: A MatrixThreadParams object is constructed
    */
    MatrixThreadParams(const MatrixSection* A, const MatrixSection* B, MatrixSection* R) {
        this->A = A;
        this->B = B;
        this->R = R;
    }
};

/*
 * description: adds two matrix sections and stores the result in a
 *                  resulting matrix section
 * return: MatrixSection*
 * precondition: 3 Matrix Sections are declared
 * postcondition: result of A + B stored in C
*/
MatrixSection* matrixAddition(const MatrixSection& A, const MatrixSection& B, MatrixSection& C) {
    for(int i = 0; i < C.dim; i++) {
        for(int j = 0; j < C.dim; j++) {
            C.M->data[i + C.rowStart][j + C.colStart] = A.M->data[i + A.rowStart][j + A.colStart] + B.M->data[i + B.rowStart][j + B.colStart];
        }
    }
    return (MatrixSection *)&C;
}

/*
 * description: subtracts two matrix sections and stores the result in a
 *                  resulting matrix section
 * return: MatrixSection*
 * precondition: 3 Matrix Sections are declared
 * postcondition: result of A - B stored in C
*/
MatrixSection* matrixSubtraction(const MatrixSection& A, const MatrixSection& B, MatrixSection& C) {
    for(int i = 0; i < C.dim; i++) {
        for(int j = 0; j < C.dim; j++) {
            C.M->data[i + C.rowStart][j + C.colStart] = A.M->data[i + A.rowStart][j + A.colStart] - B.M->data[i + B.rowStart][j + B.colStart];
        }
    }
    return (MatrixSection *)&C;
}

/*
 * description: Multiplies two matrix sections together
 * return: void *
 * precondition: A MatrixThreadParam* is declared
 * postcondition: results of A * B stored in R
*/
void* BruteForceSquareMatrixMultiplication(void* params) {
    MatrixThreadParams* matrixParams = (MatrixThreadParams*)params;

    const MatrixSection* A = matrixParams->A;
    const MatrixSection* B = matrixParams->B;
    MatrixSection* R = matrixParams->R;

    int mid = R->M->dim/2;
    for(int i = 0; i < A->dim; i++) {
        for(int j = 0; j < A->dim; j++) {
            int sum = 0;

            for(int k = 0; k < R->M->dim; k++) {
                sum += A->M->data[i + A->rowStart][k + A->colStart] * B->M->data[k + B->rowStart][j + B->colStart];
            }
            R->M->data[i + R->rowStart][j + R->colStart] = sum;
        }
    }
}


/*
 * description: Helper function for Strassen that performs Threaded Multiplication
 * return: void
 * precondition: 3 matrix sections are declared
 * postcondition: results of A * B stored in C
*/
void StrassenHelper(const MatrixSection& A, const MatrixSection& B, MatrixSection& C) {
    // uses custom input from strassen to then call threaded divide and conquer
    // Going to perform A * B and store it in C
    if(A.dim == 1 && B.dim == 1) {
        C.M->data[0][0] = A.M->data[A.rowStart][A.colStart] * B.M->data[B.rowStart][B.colStart];
    } else {
        int mid = A.dim/2;
        // split matrix into 4 sections

        MatrixThreadParams* R11Param = new MatrixThreadParams(new MatrixSection(A.M, A.rowStart, A.colStart, mid), new MatrixSection(B.M, B.rowStart, B.colStart, mid), new MatrixSection(C.M, C.rowStart, C.colStart, mid));
        MatrixThreadParams* R12Param = new MatrixThreadParams(new MatrixSection(A.M, A.rowStart, A.colStart, mid), new MatrixSection(B.M, B.rowStart, B.colStart + mid, mid), new MatrixSection(C.M, C.rowStart, C.colStart + mid, mid));
        MatrixThreadParams* R21Param = new MatrixThreadParams(new MatrixSection(A.M,  A.rowStart + mid, A.colStart, mid), new MatrixSection(B.M, B.rowStart, B.colStart, mid), new MatrixSection(C.M, C.rowStart + mid, C.colStart, mid));
        MatrixThreadParams* R22Param = new MatrixThreadParams(new MatrixSection(A.M, A.rowStart + mid, A.colStart, mid),new MatrixSection(B.M, B.rowStart, B.colStart + mid, mid), new MatrixSection(C.M, C.rowStart + mid, C.colStart + mid, mid));

        pthread_t ids[4];
        // spawn 4 threads
        pthread_create(&ids[0], NULL, &BruteForceSquareMatrixMultiplication, (void *)R11Param);
        pthread_create(&ids[1], NULL, &BruteForceSquareMatrixMultiplication, (void *)R12Param);
        pthread_create(&ids[2], NULL, &BruteForceSquareMatrixMultiplication, (void *)R21Param);
        pthread_create(&ids[3], NULL, &BruteForceSquareMatrixMultiplication, (void *)R22Param);

        // join all threads
        pthread_join(ids[0], NULL);
        pthread_join(ids[1], NULL);
        pthread_join(ids[2], NULL);
        pthread_join(ids[3], NULL);

    }
}


/*
 * description: Performs brute force matrix multiplication
 * return: SquareMatrix*
 * precondition: 2 SquareMatrix objects with equal dimensions are declared
 * postcondition: result of A * B returned
*/
SquareMatrix* BruteForce(const SquareMatrix& A, const SquareMatrix& B) {
    SquareMatrix* m = new SquareMatrix(A.dim);
    int sum;
    for(int i = 0; i < m->dim; i++) {
        for(int j = 0; j < m->dim; j++) {
            sum = 0;
            for(int k = 0; k < m->dim; k++) {
                sum += A.data[i][k] * B.data[k][j];
            }
            m->data[i][j] = sum;
        }
    }
    return m;
}

/*
 * description: Performs threaded brute force multiplication
 * return: SquareMatrix*
 * precondition: 2 SquareMatrix objects with equal dimensions are declared
 * postcondition: result of A * B returned
*/
SquareMatrix* ThreadedDivideAndConquer(const SquareMatrix& A, const SquareMatrix& B) {

    SquareMatrix* resM = new SquareMatrix(A.dim);

    if(A.dim == 1 && B.dim == 1) {
        resM->data[0][0] = A.data[0][0] * B.data[0][0];
    } else {
        int mid = A.dim / 2;

        MatrixThreadParams R11Param = MatrixThreadParams(new MatrixSection((SquareMatrix *) &A, 0, 0, mid), new MatrixSection((SquareMatrix *)&B, 0, 0, mid), new MatrixSection(resM, 0, 0, mid));
        MatrixThreadParams R12Param = MatrixThreadParams(new MatrixSection((SquareMatrix *) &A, 0, 0, mid), new MatrixSection((SquareMatrix*) &B, 0, mid,mid), new MatrixSection(resM, 0, mid, mid));
        MatrixThreadParams R21Param = MatrixThreadParams(new MatrixSection((SquareMatrix *) &A, mid, 0, mid), new MatrixSection((SquareMatrix*) &B, 0, 0, mid), new MatrixSection(resM, mid, 0, mid));
        MatrixThreadParams R22Param = MatrixThreadParams(new MatrixSection((SquareMatrix *) &A, mid, 0, mid), new MatrixSection((SquareMatrix *)&B, 0, mid, mid), new MatrixSection(resM, mid, mid, mid));

        pthread_t ids[4];

        pthread_create(&ids[0], NULL, &BruteForceSquareMatrixMultiplication, (void *)&R11Param);
        pthread_create( &ids[1], NULL, &BruteForceSquareMatrixMultiplication, (void *)&R12Param);
        pthread_create(&ids[2], NULL, &BruteForceSquareMatrixMultiplication, (void *)&R21Param);
        pthread_create(&ids[3], NULL, &BruteForceSquareMatrixMultiplication, (void *)&R22Param);

        pthread_join(ids[0], NULL);
        pthread_join(ids[1], NULL);
        pthread_join(ids[2], NULL);
        pthread_join(ids[3], NULL);
    }
    return resM;

}

/*
 * description: Performs Strassen matrix multiplication on two matrices
 * return: SquareMatrix*
 * precondition: 2 SquareMatrix objects with equal dimensions are declared
 * postcondition: result of A * B returned
*/
SquareMatrix* Strassen(const SquareMatrix& A, const SquareMatrix& B) {

    int mid = A.dim/2;

    // no new matrices being built, just sections being assigned
    MatrixSection A11 = MatrixSection((SquareMatrix *)&A, 0, 0, mid);
    MatrixSection A12 = MatrixSection((SquareMatrix *)&A, 0, mid, mid);
    MatrixSection A21 = MatrixSection((SquareMatrix *)&A, mid, 0, mid);
    MatrixSection A22 = MatrixSection((SquareMatrix *)&A, mid, mid, mid);
    MatrixSection B11 = MatrixSection((SquareMatrix *)&B, 0, 0, mid);
    MatrixSection B12 = MatrixSection((SquareMatrix *)&B, 0, mid, mid);
    MatrixSection B21 = MatrixSection((SquareMatrix *)&B, mid, 0, mid);
    MatrixSection B22 = MatrixSection((SquareMatrix *)&B, mid, mid, mid);

    MatrixSection* sMatrices[7];

    // only time that this function will change with input size. No way around it
    for(int i = 0; i < 7; i++) {
        sMatrices[i] = new MatrixSection(new SquareMatrix(mid), 0, 0, mid); // death
    }

    SquareMatrix* resM = new SquareMatrix(A.dim);

    // storing results of operations
    MatrixSection storeRes1 = MatrixSection(new SquareMatrix(mid), 0, 0, mid);
    MatrixSection storeRes2 = MatrixSection(new SquareMatrix(mid), 0, 0, mid);

    // s1
    StrassenHelper(A11, *matrixSubtraction(B12, B22, storeRes1), *sMatrices[0]);
    // s2
    StrassenHelper(*matrixAddition(A11, A12, storeRes1), B22, *sMatrices[1]);
    // s3
    StrassenHelper(*matrixAddition(A21, A22, storeRes1), B11, *sMatrices[2]);
    // s4
    StrassenHelper(A22, *matrixSubtraction(B21, B11, storeRes1), *sMatrices[3]);
    // s5
    StrassenHelper(*matrixAddition(A11, A22, storeRes1), *matrixAddition(B11, B22, storeRes2), *sMatrices[4]);
    // s6
    StrassenHelper(*matrixSubtraction(A12, A22, storeRes1), *matrixAddition(B21, B22, storeRes2), *sMatrices[5]);
    // s7
    StrassenHelper(*matrixSubtraction(A11, A21, storeRes1), *matrixAddition(B11, B12, storeRes2), *sMatrices[6]);


    MatrixSection I = MatrixSection(resM, 0, 0, mid);
    MatrixSection J = MatrixSection(resM, 0, mid, mid);
    MatrixSection K = MatrixSection(resM, mid, 0, mid);
    MatrixSection L = MatrixSection(resM, mid, mid, mid);

    // I
    matrixAddition(*sMatrices[4], *sMatrices[3], storeRes1);
    matrixSubtraction(storeRes1, *sMatrices[1], storeRes2);
    matrixAddition(storeRes2, *sMatrices[5], I);

    // J
    matrixAddition(*sMatrices[0], *sMatrices[1], J);

    // K
    matrixAddition(*sMatrices[2], *sMatrices[3], K);

    // L
    matrixAddition(*sMatrices[0], *sMatrices[4], storeRes1);
    matrixSubtraction(storeRes1, *sMatrices[2], storeRes2);
    matrixSubtraction(storeRes2, *sMatrices[6], L);
    return resM;

}

/*
 * description: Performs threaded Strassen matrix multiplication on two matrices
 * return: SquareMatrix*
 * precondition: 2 SquareMatrix objects with equal dimensions are declared
 * postcondition: result of A * B returned
*/
SquareMatrix* ThreadedStrassen(const SquareMatrix& A, const SquareMatrix& B) {

    int mid = A.dim/2;

    // no new matrices being built, just sections being assigned
    MatrixSection A11 =  MatrixSection((SquareMatrix *)&A, 0, 0, mid);
    MatrixSection A12 = MatrixSection((SquareMatrix *)&A, 0, mid, mid);
    MatrixSection A21 = MatrixSection((SquareMatrix *)&A, mid, 0, mid);
    MatrixSection A22 = MatrixSection((SquareMatrix *)&A, mid, mid, mid);
    MatrixSection B11 = MatrixSection((SquareMatrix *)&B, 0, 0, mid);
    MatrixSection B12 = MatrixSection((SquareMatrix *)&B, 0, mid, mid);
    MatrixSection B21 = MatrixSection((SquareMatrix *)&B, mid, 0, mid);
    MatrixSection B22 = MatrixSection((SquareMatrix *)&B, mid, mid, mid);

    MatrixSection* sMatrices[7];

    // only time that this function will change with input size. No way around it
    for(int i = 0; i < 7; i++) {
        sMatrices[i] = new MatrixSection(new SquareMatrix(mid), 0, 0, mid); // death
    }

    SquareMatrix* resM = new SquareMatrix(A.dim);

    // storing results of operations
    MatrixSection storeRes1 = MatrixSection(new SquareMatrix(mid), 0, 0, mid);
    MatrixSection storeRes2 = MatrixSection(new SquareMatrix(mid), 0, 0, mid);
    MatrixSection storeRes3 = MatrixSection(new SquareMatrix(mid), 0, 0, mid);
    MatrixSection storeRes4 = MatrixSection(new SquareMatrix(mid), 0, 0, mid);
    MatrixSection storeRes5 = MatrixSection(new SquareMatrix(mid), 0, 0, mid);
    MatrixSection storeRes6 = MatrixSection(new SquareMatrix(mid), 0, 0, mid);


    pthread_t ids[4];
    MatrixThreadParams s1 = MatrixThreadParams(&A11, matrixSubtraction(B12, B22, storeRes1), sMatrices[0]);
    MatrixThreadParams s2 = MatrixThreadParams(matrixAddition(A11, A12, storeRes2), &B22, sMatrices[1]);
    MatrixThreadParams s3 = MatrixThreadParams(matrixAddition(A21, A22, storeRes3), &B11, sMatrices[2]);
    MatrixThreadParams s4 = MatrixThreadParams(&A22, matrixSubtraction(B21, B11, storeRes4), sMatrices[3]);

    pthread_create(&ids[0] , NULL, &BruteForceSquareMatrixMultiplication, (void *)&s1);
    pthread_create(&ids[1] , NULL, &BruteForceSquareMatrixMultiplication, (void *)&s2);
    pthread_create(&ids[2] , NULL, &BruteForceSquareMatrixMultiplication, (void *)&s3);
    pthread_create(&ids[3] , NULL, &BruteForceSquareMatrixMultiplication, (void *)&s4);

    pthread_join(ids[0], NULL);
    pthread_join(ids[1], NULL);
    pthread_join(ids[2], NULL);
    pthread_join(ids[3], NULL);

    MatrixThreadParams s5 = MatrixThreadParams(matrixAddition(A11, A22, storeRes1), matrixAddition(B11, B22, storeRes2), sMatrices[4]);
    MatrixThreadParams s6 = MatrixThreadParams(matrixSubtraction(A12, A22, storeRes3), matrixAddition(B21, B22, storeRes4), sMatrices[5]);
    MatrixThreadParams s7 = MatrixThreadParams(matrixSubtraction(A11, A21, storeRes5), matrixAddition(B11, B12, storeRes6), sMatrices[6]);

    pthread_create(&ids[0] , NULL, &BruteForceSquareMatrixMultiplication, (void *)&s5);
    pthread_create(&ids[1] , NULL, &BruteForceSquareMatrixMultiplication, (void *)&s6);
    pthread_create(&ids[2] , NULL, &BruteForceSquareMatrixMultiplication, (void *)&s7);

    pthread_join(ids[0], NULL);
    pthread_join(ids[1], NULL);
    pthread_join(ids[2], NULL);

    MatrixSection *I, *J, *K, *L;
    I = new MatrixSection(resM, 0, 0, mid);
    J = new MatrixSection(resM, 0, mid, mid);
    K = new MatrixSection(resM, mid, 0, mid);
    L = new MatrixSection(resM, mid, mid, mid);

    // I
    matrixAddition(*sMatrices[4], *sMatrices[3], storeRes1);
    matrixSubtraction(storeRes1, *sMatrices[1], storeRes2);
    matrixAddition(storeRes2, *sMatrices[5], *I);

    // J
    matrixAddition(*sMatrices[0], *sMatrices[1], *J);

    // K
    matrixAddition(*sMatrices[2], *sMatrices[3], *K);

    // L
    matrixAddition(*sMatrices[0], *sMatrices[4], storeRes1);
    matrixSubtraction(storeRes1, *sMatrices[2], storeRes2);
    matrixSubtraction(storeRes2, *sMatrices[6], *L);
    return resM;
}

#endif//MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H
