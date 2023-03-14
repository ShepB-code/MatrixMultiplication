#include <iostream>
#include <thread>
#include "ThreadedMatrixMult.h"
using namespace std;


ostream& printMatrix(const SquareMatrix&, ostream&);
SquareMatrix* readMatrix(istream&);

int main() {


    //read in matrices
    SquareMatrix* m1 = readMatrix(cin);
    SquareMatrix* m2 = readMatrix(cin);

    cout << "---------DATA---------" << endl;
    cout << "Matrix 1" << endl;
    printMatrix(*m1, cout);
    cout << "Matrix 2" << endl;
    printMatrix(*m2, cout);
    cout << "------------------" << endl;

    cout << "---------BRUTE FORCE---------" << endl;
    SquareMatrix* bruteForce = BruteForce(*m1, *m2);
    printMatrix(*bruteForce, cout);

    cout << "---------D&CThreaded---------" << endl;
    SquareMatrix* dCThread = ThreadedDivideAndConquer(*m1, *m2);
    printMatrix(*dCThread, cout);

    cout << "---------Strassen---------" << endl;
    SquareMatrix* strass = Strassen(*m1, *m2);
    printMatrix(*strass, cout);
    cout << "---------ThreadedStrassen---------" << endl;


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
    for(int i = 0; i < m->dim; i++) {
        for(int j = 0; j < m->dim; j++) {
            in >> m->data[i][j];
        }
    }
    return m;
}



