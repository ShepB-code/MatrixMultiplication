#include "ThreadedMatrixMult.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>
using namespace std;
using namespace chrono;

const int DIMENSION = 1024;
const int tests = 20;

ostream& printMatrix(const SquareMatrix&, ostream&);
SquareMatrix* readMatrix(istream&);
void makeMatrix(SquareMatrix&);

int main() {

    SquareMatrix testA(DIMENSION), testB(DIMENSION);

    makeMatrix(testA);
    makeMatrix(testB);
    {
        cout << "---------BRUTE FORCE---------" << endl;
        vector<decltype(high_resolution_clock::now() - high_resolution_clock::now())> data;
        for (int m = 0; m < tests; m++) {
            auto start = high_resolution_clock::now();

            BruteForce(testA, testB);

            auto stop = high_resolution_clock::now();
            data.emplace_back(duration_cast<milliseconds>(stop - start));
        }
        double average = 0;

        for (auto i: data) {
            average += i.count();
        }
        average /= data.size();

        cout << setprecision(20) << average << endl;
    }

    {
        cout << "---------D&CThreaded---------" << endl;
        vector<decltype(high_resolution_clock::now() - high_resolution_clock::now())> data;
        for (int m = 0; m < tests; m++) {
            auto start = high_resolution_clock::now();

            ThreadedDivideAndConquer(testA, testB);

            auto stop = high_resolution_clock::now();
            data.emplace_back(duration_cast<milliseconds>(stop - start));
        }
        double average = 0;

        for (auto i: data) {
            average += i.count();
        }
        average /= data.size();

        cout << setprecision(20) << average << endl;
    }

    {
        cout << "---------Strassen---------" << endl;
        vector<decltype(high_resolution_clock::now() - high_resolution_clock::now())> data;
        for (int m = 0; m < tests; m++) {
            auto start = high_resolution_clock::now();

            Strassen(testA, testB);

            auto stop = high_resolution_clock::now();
            data.emplace_back(duration_cast<milliseconds>(stop - start));
        }
        double average = 0;

        for (auto i: data) {
            average += i.count();
        }
        average /= data.size();

        cout << setprecision(20) << average << endl;
    }

    {
        cout << "---------ThreadedStrassen---------" << endl;

        vector<decltype(high_resolution_clock::now() - high_resolution_clock::now())> data;
        for (int m = 0; m < tests; m++) {
            auto start = high_resolution_clock::now();

            ThreadedStrassen(testA, testB);

            auto stop = high_resolution_clock::now();
            data.emplace_back(duration_cast<milliseconds>(stop - start));
        }
        double average = 0;

        for (auto i: data) {
            average += i.count();
        }
        average /= data.size();

        cout << setprecision(20) << average << endl;
    }

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

void makeMatrix(SquareMatrix& m) {
    for (int i = 0; i < m.dim; i++) {
        for (int j = 0; j < m.dim; j++) {
            m.data[i][j] = rand() % 100;
        }
    }
}


