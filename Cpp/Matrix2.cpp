#include "../Headers/Matrices.h"
#include <iostream>
#include <fstream>
using namespace std;
int main(){
    int m, n;
    cout << " Enter size (rows x columns): ";
    cin >> m >> n;
    Matrix<float> M(m, n);
    cin >> M;
    cout << " The matrix = \n";
    cout << M << endl;
    Matrix<float> M12 = (M * M);
    for(int i = 0; i < 4; i++){
        cout << " M^"<<i+2<<" = \n";
        cout << M12 <<endl;
        M12 = M12 * M;
    }
}