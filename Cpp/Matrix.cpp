#include "../Headers/Matrices.h"
#include <iostream>
#include <fstream>
using namespace std;
int main(){
    Matrix<float> M(4, 4);
    ifstream ifs("Matrix2.dat");
    if(!ifs)
    {
        ofstream ofs("Matrix2.dat");
        cin >> M;
        ofs << M;
        cout << M;
    }
    else 
    {
        cout << " reading from file\n";
        ifs >> M;
        cout << M;
        cout << " determinant of matrix = "<<M.determinant() << endl;
        cout << " cofactor of matrix = \n" <<M.cofactor() << endl;
        cout << " adjacent of matrix = \n" <<M.adjacent() << endl;
    }
}