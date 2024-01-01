#include "../Headers/Matrices.h"
#include <iostream>
#include <fstream>
using namespace std;
int main(){
    int m1, n1, m2, n2;
    cout << " Enter size 1 (rows x columns): ";
    cin >> m1 >> n1;
    cout << " Enter size 2 (rows x columns): ";
    cin >> m2 >> n2;
    Matrix<float> M1(m1, n1);
    Matrix<float> M2(m2, n2);
    cout << " Enter matrix 1 entries: \n";
    cin >> M1;
    cout << " Enter matrix 2 entries: \n";
    cin >> M2;
    cout << " Matrix 1 = \n";
    cout << M1 << endl;
    cout << " Matrix 2 = \n";
    cout << M2 << endl;
    cout << " Matrix 1 * Matrix 2 = \n";
    cout << M1 * M2 << endl;
}