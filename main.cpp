#include <iostream>
#include <thread>
#include <assert.h>
#include <chrono>
#include "ThreadedMatrixMult.h"
using namespace std;

const int DIMENSION = 1024;

ostream& printMatrix(const SquareMatrix&, ostream&);
SquareMatrix* readMatrix(istream&);
SquareMatrix* makeMatrix(int);
bool compareMatrices(const SquareMatrix&, const SquareMatrix&);

int main() {
    //read in matrices through file
    /*
    cout << "FIRST MATRIX" << endl;
    SquareMatrix* m1 = readMatrix(cin);
    cout << "SECOND MATRIX" << endl;
    SquareMatrix* m2 = readMatrix(cin);
     */
    SquareMatrix* m1 = makeMatrix(DIMENSION);
    SquareMatrix* m2 = makeMatrix(DIMENSION);
    //cout << "---------DATA---------" << endl;
    //cout << "Matrix 1" << endl;
    //printMatrix(*m1, cout);
    //cout << "Matrix 2" << endl;
    //printMatrix(*m2, cout);


    cout << "---------BRUTE FORCE---------" << endl;
    auto start_brute = chrono::high_resolution_clock::now();
    SquareMatrix* bruteForce = BruteForce(*m1, *m2);
    auto stop_brute = chrono::high_resolution_clock::now();
    //printMatrix(*bruteForce, cout);

    cout << "---------D&CThreaded---------" << endl;

    auto start_dCThread = chrono::high_resolution_clock::now();
    SquareMatrix* dCThread = ThreadedDivideAndConquer(*m1, *m2);
    auto stop_dCThread = chrono::high_resolution_clock::now();
    //printMatrix(*dCThread, cout);


    cout << "---------Strassen---------" << endl;
    auto start_strass = chrono::high_resolution_clock::now();
    SquareMatrix* strass = Strassen(*m1, *m2);
    auto stop_strass = chrono::high_resolution_clock::now();
    //printMatrix(*strass, cout);

    cout << "---------ThreadedStrassen---------" << endl;
    auto start_tStrass = chrono::high_resolution_clock::now();
    SquareMatrix* tStrass = ThreadedStrassen(*m1, *m2);
    auto stop_tStrass = chrono::high_resolution_clock::now();
    //printMatrix(*tStrass, cout);


    cout << "---------COMPARING---------" << endl;
    if(compareMatrices(*bruteForce, *dCThread) && compareMatrices(*bruteForce, *strass) &&
            compareMatrices(*bruteForce, *tStrass) ) { //&& compareMatrices(*bruteForce, *strass)) {
        cout << "All matrices are the same" << endl;
    }


    cout << "---------" << DIMENSION << "x" << DIMENSION << " EXECUTION TIME---------" << endl;
    cout << "Brute Force: " << chrono::duration_cast<chrono::microseconds>(stop_brute- start_brute).count() << endl;
    cout << "dCThread: " << chrono::duration_cast<chrono::microseconds>(stop_dCThread - start_dCThread).count() << endl;
    cout << "Strassen: " << chrono::duration_cast<chrono::microseconds>(stop_strass- start_strass).count() << endl;
    cout << "Threaded Strassen: " << chrono::duration_cast<chrono::microseconds>(stop_tStrass- start_tStrass).count() << endl;

    //  SquareMatrix* divideAndConquer = ThreadedDivideAndConquer(m1, m2);




    // perform brute force

    return 0;
}

ostream& printMatrix(const SquareMatrix& m, ostream& out) {
    for(int i = 0; i < m.dim; i++) {
        for(int j = 0; j < m.dim; j++) {
            out << m.data[i][j] << " ";
        }
        out << endl;
    }
    return out;
}

SquareMatrix* readMatrix(istream& in) {
    int dim;
    cin >> dim;
    SquareMatrix* m = new SquareMatrix(dim);
    int count = 0;
    for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
            in >> m->data[i][j];
            cout << "Curr: " << m->data[i][j] << " " << "Count: " << ++count << endl;
        }
    }
    return m;
}

SquareMatrix* makeMatrix(int dim) {
    SquareMatrix* m = new SquareMatrix(dim);
    int count = 1;
    for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
            //m->data[i][j] = rand() % 100;
            m->data[i][j] = count;
        }
        count++;
    }
    return m;
}

bool compareMatrices(const SquareMatrix& A, const SquareMatrix& B) {
    bool same = true;
    for(int i = 0; i < A.dim && same; i++) {
        for(int j = 0; j < A.dim && same; j++) {
            same = A.data[i][j] == B.data[i][j];
        }
    }
    assert(same);
    return same;
}

