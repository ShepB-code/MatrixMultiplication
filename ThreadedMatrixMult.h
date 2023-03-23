//
// Created by shepa on 3/10/2023.
//

#ifndef MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H
#define MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H

#include <pthread.h>
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
    ~SquareMatrix() {
        for(int i = 0; i < dim; i++) {
            delete[] data;
        }
    }
};

struct MatrixSection {
    // refers to a complete matrix, and knows its bounds
    SquareMatrix* M;
    int rowStart, colStart, dim;
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
    bool all;
    MatrixThreadParams(const MatrixSection* A, const MatrixSection* B, MatrixSection* R, bool all) {
        this->A = A;
        this->B = B;
        this->R = R;
        this->all = all;
    }
    //MatrixThreadParams(const SquareMatrix *a, const SquareMatrix *b, SquareMatrix *r, int rowEnd, int colEnd, int dim, int rowStart, int colStart) : A(a), B(b), R(r), rowEnd(rowEnd), colEnd(colEnd), dim(dim), rowStart(rowStart), colStart(colStart) {}
};

MatrixSection* matrixAddition(const MatrixSection& A, const MatrixSection& B, MatrixSection& C) {
    //SquareMatrix* C = new SquareMatrix(A.dim);
    //cout << A.dim << "x" << B.dim << "->" << C.dim << endl;

    for(int i = 0; i < C.dim; i++) {
        for(int j = 0; j < C.dim; j++) {
            C.M->data[i + C.rowStart][j + C.colStart] = A.M->data[i + A.rowStart][j + A.colStart] + B.M->data[i + B.rowStart][j + B.colStart];
        }
    }
    return (MatrixSection *)&C;
}

MatrixSection* matrixSubtraction(const MatrixSection& A, const MatrixSection& B, MatrixSection& C) {
    //cout << A.dim << "x" << B.dim << "->" << C.dim << endl;
    for(int i = 0; i < C.dim; i++) {
        for(int j = 0; j < C.dim; j++) {
            C.M->data[i + C.rowStart][j + C.colStart] = A.M->data[i + A.rowStart][j + A.colStart] - B.M->data[i + B.rowStart][j + B.colStart];
        }
    }
    return (MatrixSection *)&C;
}

void* BruteForceSquareMatrixMultiplication(void* params) {
    MatrixThreadParams* matrixParams = (MatrixThreadParams*)params;

    //SquareMatrix* resM = new SquareMatrix{matrixPackage->A->dim, new int*[matrixPackage->A->dim]};

    const MatrixSection* A = matrixParams->A;
    const MatrixSection* B = matrixParams->B;
    MatrixSection* R = matrixParams->R;

    //cout << "A (" << A->rowStart << "," << A->colStart << ") " << " B ("<< B->rowStart << "," << B->colStart << ") " << "R (" << R->rowStart << "," << R->colStart << ")" << endl;


    int mid = R->M->dim/2;
    for(int i = 0; i < mid; i++) {
        for(int j = 0; j < mid; j++) {
            //matrixParams->R->data[i][j] = 0;
            int sum = 0;

            int kBound = matrixParams->all ? R->M->dim : min(A->M->dim, B->M->dim);
            for(int k = 0; k < kBound; k++) {
                if(matrixParams->all) {
                    sum += A->M->data[i][k] * B->M->data[k][j];
                } else {
                    if (A->M->dim > B->M->dim) {
                        sum += A->M->data[i + A->rowStart][k + A->colStart] * B->M->data[k][j];
                    } else if (B->M->dim > A->M->dim) {
                        sum += A->M->data[i][k] * B->M->data[k + B->rowStart][j + B->colStart];
                    } else {
                        sum += A->M->data[i][k] * B->M->data[k][j];
                    }
                }
            }
            R->M->data[i + R->rowStart][j + R->colStart] = sum;
        }
    }
    //cout << "Done" << endl;
    //pthread_exit(NULL);
    //return params;
}

void StrassenHelper(const MatrixSection& A, const MatrixSection& B, MatrixSection& C) {
    // uses custom input from strassen to then call threaded divide and conquer
    // Going to perform A * B and store it in C
    if(A.dim == 1 && B.dim == 1) {
        C.M->data[0][0] = A.M->data[A.rowStart][A.colStart] * B.M->data[B.rowStart][B.colStart];
    } else {
        int mid = A.dim/2;
        // split matrix into 4 sections
        /*
        MatrixThreadParams* R11Param = new MatrixThreadParams(new MatrixSection(A.M, A.rowStart, A.colStart, mid), new MatrixSection(B.M, B.rowStart, B.colStart, mid), new MatrixSection(C.M, C.rowStart, C.colStart, mid), false);
        MatrixThreadParams* R12Param = new MatrixThreadParams(new MatrixSection(A.M, A.rowStart, A.colStart + mid, mid), new MatrixSection(B.M, B.rowStart, B.colStart + mid, mid), new MatrixSection(C.M, C.rowStart, C.colStart + mid, mid), false);
        MatrixThreadParams* R21Param = new MatrixThreadParams(new MatrixSection(A.M,  A.rowStart + mid, A.colStart, mid), new MatrixSection(B.M, B.rowStart + mid, B.colStart, mid), new MatrixSection(C.M, C.rowStart + mid, C.colStart, mid), false);
        MatrixThreadParams* R22Param = new MatrixThreadParams(new MatrixSection(A.M, A.rowStart + mid, A.colStart + mid, mid), new MatrixSection(B.M, B.rowStart + mid, B.colStart + mid, mid), new MatrixSection(C.M, C.rowStart + mid, C.colStart + mid, mid), false);

        */
        MatrixThreadParams* R11Param = new MatrixThreadParams(new MatrixSection(A.M, 0, 0, mid), new MatrixSection(B.M, 0, 0, mid), new MatrixSection(C.M, C.rowStart, C.colStart, mid), false);
        MatrixThreadParams* R12Param = new MatrixThreadParams(new MatrixSection(A.M, 0,  mid, mid), new MatrixSection(B.M, 0, mid, mid), new MatrixSection(C.M, C.rowStart, C.colStart + mid, mid), false);
        MatrixThreadParams* R21Param = new MatrixThreadParams(new MatrixSection(A.M,  mid, 0, mid), new MatrixSection(B.M,  mid, 0, mid), new MatrixSection(C.M, C.rowStart + mid, C.colStart, mid), false);
        MatrixThreadParams* R22Param = new MatrixThreadParams(new MatrixSection(A.M, mid, mid, mid), new MatrixSection(B.M, mid, mid, mid), new MatrixSection(C.M, C.rowStart + mid, C.colStart + mid, mid), false);

        pthread_t id1, id2, id3, id4;
        // spawn 4 threads


        BruteForceSquareMatrixMultiplication((void*)R11Param);
        BruteForceSquareMatrixMultiplication((void*)R12Param);
        BruteForceSquareMatrixMultiplication((void*)R21Param);
        BruteForceSquareMatrixMultiplication((void*)R22Param);

        /*
        pthread_create(&id1, NULL, &BruteForceSquareMatrixMultiplication, (void *)R11Param);
        pthread_create(&id2, NULL, &BruteForceSquareMatrixMultiplication, (void *)R12Param);
        pthread_create(&id3, NULL, &BruteForceSquareMatrixMultiplication, (void *)R21Param);
        pthread_create(&id4, NULL, &BruteForceSquareMatrixMultiplication, (void *)R22Param);


        // join all threads
        pthread_join(id1, NULL);
        pthread_join(id2, NULL);
        pthread_join(id3, NULL);
        pthread_join(id4, NULL);
        */
    }
}


// BOOTH FUNCTIONS (defined by me)
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
SquareMatrix* ThreadedDivideAndConquer(const SquareMatrix& A, const SquareMatrix& B) {

    SquareMatrix* resM = new SquareMatrix(A.dim);

    if(A.dim == 1 && B.dim == 1) {
        resM->data[0][0] = A.data[0][0] * B.data[0][0];
    } else {
        int mid = A.dim / 2;

        MatrixThreadParams* R11Param = new MatrixThreadParams(new MatrixSection((SquareMatrix *) &A, 0, 0, mid), new MatrixSection((SquareMatrix *)&B, 0, 0, mid), new MatrixSection(resM, 0, 0, mid), true);
        MatrixThreadParams* R12Param = new MatrixThreadParams(new MatrixSection((SquareMatrix *) &A, 0, mid, mid), new MatrixSection((SquareMatrix*) &B, 0, mid,mid), new MatrixSection(resM, 0, mid, mid), true);
        MatrixThreadParams* R21Param = new MatrixThreadParams(new MatrixSection((SquareMatrix *) &A, mid, 0, mid), new MatrixSection((SquareMatrix*) &B, mid, 0, mid), new MatrixSection(resM, mid, 0, mid), true);
        MatrixThreadParams* R22Param = new MatrixThreadParams(new MatrixSection((SquareMatrix *) &A, mid, mid, mid), new MatrixSection((SquareMatrix *)&B, mid, mid, mid), new MatrixSection(resM, mid, mid, mid),true);

        pthread_t id1, id2, id3, id4;


        cout << "R11: A(" << R11Param->A->rowStart << "," << R11Param->A->colStart << ") B(" << R11Param->B->rowStart << "," << R11Param->B->colStart << ") R(" << R11Param->R->rowStart << "," << R11Param->R->colStart << ")" << endl;
        cout << "R12: A(" << R12Param->A->rowStart << "," << R12Param->A->colStart << ") B(" << R12Param->B->rowStart << "," << R12Param->B->colStart << ") R(" << R12Param->R->rowStart << "," << R12Param->R->colStart << ")" << endl;
        cout << "R21: A(" << R21Param->A->rowStart << "," << R21Param->A->colStart << ") B(" << R21Param->B->rowStart << "," << R21Param->B->colStart << ") R(" << R21Param->R->rowStart << "," << R21Param->R->colStart << ")" << endl;
        cout << "R22: A(" << R22Param->A->rowStart << "," << R22Param->A->colStart << ") B(" << R22Param->B->rowStart << "," << R22Param->B->colStart << ") R(" << R22Param->R->rowStart << "," << R22Param->R->colStart << ")" << endl;

        /*
        BruteForceSquareMatrixMultiplication((void *)R11Param);
        BruteForceSquareMatrixMultiplication((void *)R12Param);
        BruteForceSquareMatrixMultiplication((void *)R21Param);
        BruteForceSquareMatrixMultiplication((void *)R22Param);
        */

        pthread_create(&id1, NULL, &BruteForceSquareMatrixMultiplication, (void *)R11Param);
        pthread_create( &id2, NULL, &BruteForceSquareMatrixMultiplication, (void *)R12Param);
        pthread_create(&id3, NULL, &BruteForceSquareMatrixMultiplication, (void *)R21Param);
        pthread_create(&id4, NULL, &BruteForceSquareMatrixMultiplication, (void *)R22Param);

        pthread_join(id1, NULL);
        pthread_join(id2, NULL);
        pthread_join(id3, NULL);
        pthread_join(id4, NULL);
    }
    return resM;

}
SquareMatrix* Strassen(const SquareMatrix& A, const SquareMatrix& B) {

    int mid = A.dim/2;

    // no new matrices being built, just sections being assigned
    MatrixSection* A11 = new MatrixSection((SquareMatrix *)&A, 0, 0, mid);
    MatrixSection* A12 = new MatrixSection((SquareMatrix *)&A, 0, mid, mid);
    MatrixSection* A21 = new MatrixSection((SquareMatrix *)&A, mid, 0, mid);
    MatrixSection* A22 = new MatrixSection((SquareMatrix *)&A, mid, mid, mid);
    MatrixSection* B11 = new MatrixSection((SquareMatrix *)&B, 0, 0, mid);
    MatrixSection* B12 = new MatrixSection((SquareMatrix *)&B, 0, mid, mid);
    MatrixSection* B21 = new MatrixSection((SquareMatrix *)&B, mid, 0, mid);
    MatrixSection* B22 = new MatrixSection((SquareMatrix *)&B, mid, mid, mid);

    vector<MatrixSection *> sMatrices(7);

    // only time that this function will change with input size. No way around it
    for(int i = 0; i < sMatrices.size(); i++) {
        sMatrices.at(i) = new MatrixSection(new SquareMatrix(mid), 0, 0, mid); // death
    }

    SquareMatrix* resM = new SquareMatrix(A.dim);

    // storing results of operations
    MatrixSection* storeRes1 = new MatrixSection(new SquareMatrix(mid), 0, 0, mid);
    MatrixSection* storeRes2 = new MatrixSection(new SquareMatrix(mid), 0, 0, mid);

    // s1
    StrassenHelper(*A11, *matrixSubtraction(*B12, *B22, *storeRes1), *sMatrices.at(0));
    // s2
    StrassenHelper(*matrixAddition(*A11, *A12, *storeRes1), *B22, *sMatrices.at(1));
    // s3
    StrassenHelper(*matrixAddition(*A21, *A22, *storeRes1), *B11, *sMatrices.at(2));
    // s4
    StrassenHelper(*A22, *matrixSubtraction(*B21, *B11, *storeRes1), *sMatrices.at(3));
    // s5
    StrassenHelper(*matrixAddition(*A11, *A22, *storeRes1), *matrixAddition(*B11, *B22, *storeRes2), *sMatrices.at(4));
    // s6
    StrassenHelper(*matrixSubtraction(*A12, *A22, *storeRes1), *matrixAddition(*B21, *B22, *storeRes2), *sMatrices.at(5));
    // s7
    StrassenHelper(*matrixSubtraction(*A11, *A21, *storeRes1), *matrixAddition(*B11, *B12, *storeRes2), *sMatrices.at(6));

    MatrixSection *I, *J, *K, *L;
    I = new MatrixSection(resM, 0, 0, mid);
    J = new MatrixSection(resM, 0, mid, mid);
    K = new MatrixSection(resM, mid, 0, mid);
    L = new MatrixSection(resM, mid, mid, mid);

    cout << "BIG BODY L" << endl;


    // I
    matrixAddition(*sMatrices.at(4 ), *sMatrices.at(3), *storeRes1);
    matrixSubtraction(*storeRes1, *sMatrices.at(1), *storeRes2);
    matrixAddition(*storeRes2, *sMatrices.at(5), *I);
    //matrixSubtraction(*sMatrices.at(3), *sMatrices.at(1), *storeRes1);
    //matrixAddition(*sMatrices.at(5), *storeRes1, *storeRes2);
    //matrixAddition(*sMatrices.at(4), *storeRes2, *I);

    // J
    matrixAddition(*sMatrices.at(0), *sMatrices.at(1), *J);

    // K
    matrixAddition(*sMatrices.at(2), *sMatrices.at(3), *K);

    // L
    matrixAddition(*sMatrices.at(0), *sMatrices.at(4), *storeRes1);
    matrixSubtraction(*storeRes1, *sMatrices.at(2), *storeRes2);
    matrixSubtraction(*storeRes2, *sMatrices.at(6), *L);
    //MatrixSection* helpMe = matrixAddition(*sMatrices.at(4),  *matrixAddition(*sMatrices.at(5), , *I);
    /*
    SquareMatrix* I = matrixAddition(*s5, *matrixAddition(*s6, *matrixSubtraction(*s4, *s2))); // (s5+(s6+(s4-s2)))
    SquareMatrix* J = matrixAddition(*s1, *s2); // (s1 + s2)
    SquareMatrix* K = matrixAddition(*s3, *s4); // (s3 + s4)
    SquareMatrix* L = matrixAddition(*matrixSubtraction(*matrixSubtraction(*s1, *s7), *s3), *s5); // (s1-(s7-(s3+s5)))
    */
    return resM;

}
SquareMatrix* ThreadedStrassen(const SquareMatrix& A, const SquareMatrix& B) {

    int mid = A.dim/2;
    SquareMatrix* A11 = new SquareMatrix(mid); // A
    SquareMatrix* A12 = new SquareMatrix(mid); // B
    SquareMatrix* A21 = new SquareMatrix(mid); // C
    SquareMatrix* A22 = new SquareMatrix(mid); // D
    SquareMatrix* B11 = new SquareMatrix(mid); // E
    SquareMatrix* B12 = new SquareMatrix(mid); // F
    SquareMatrix* B21 = new SquareMatrix(mid); // G
    SquareMatrix* B22 = new SquareMatrix(mid); // H

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

    /*
    SquareMatrix* s1 = ThreadedDivideAndConquer(*A11, *matrixSubtraction(*B12, *B22)); // A(F-H)
    SquareMatrix* s2 = ThreadedDivideAndConquer(*matrixAddition(*A11, *A12), *B22); // (A+B)F
    SquareMatrix* s3 = ThreadedDivideAndConquer(*matrixAddition(*A21, *A22), *B11); // (C+D)E
    SquareMatrix* s4 = ThreadedDivideAndConquer(*A22, *matrixSubtraction(*B21, *B11)); // D(G-E)
    SquareMatrix* s5 = ThreadedDivideAndConquer(*matrixAddition(*A11, *A22), *matrixAddition(*B11, *B22)); // (A+D)(E+H)
    SquareMatrix* s6 = ThreadedDivideAndConquer(*matrixSubtraction(*A12, *A22), *matrixAddition(*B21, *B22)); // (B-D)(G+H)
    SquareMatrix* s7 = ThreadedDivideAndConquer(*matrixSubtraction(*A11, *A21), *matrixAddition(*B11, *B12)); // (A-C)(E+F)

    SquareMatrix* I = matrixAddition(*s5, *matrixAddition(*s6, *matrixSubtraction(*s4, *s2))); // (s5+(s6+(s4-s2)))
    SquareMatrix* J = matrixAddition(*s1, *s2); // (s1 + s2)
    SquareMatrix* K = matrixAddition(*s3, *s4); // (s3 + s4)
    SquareMatrix* L = matrixAddition(*matrixSubtraction(*matrixSubtraction(*s1, *s7), *s3), *s5); // (s1-(s7-(s3+s5)))

    SquareMatrix* resM = new SquareMatrix(A.dim);

    for(int i = 0; i < mid; i++) {
        for(int j = 0; j < mid; j++) {
            resM->data[i][j] = I->data[i][j];
            resM->data[i][j + mid] = J->data[i][j];
            resM->data[i + mid][j] = K->data[i][j];
            resM->data[i + mid][j + mid] = L->data[i][j];
        }
    }
     */
     SquareMatrix* resM = new SquareMatrix(A.dim);

    return resM;
}




#endif//MATRIXMULTIPLICATION_THREADEDMATRIXMULT_H
