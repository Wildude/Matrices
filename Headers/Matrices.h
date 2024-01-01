#include <iostream>
#include <fstream>
#include <limits>
using namespace std;
typedef unsigned short ushort;
template <class T>
class Matrix{
    T** M;
    ushort m, n;
    void make_copy(Matrix mat){
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++)
            M[i][j] = mat.M[i][j];
        }
    }
    public:
    /*
    }
    */
    T& giveEntry(int ithRow, int jthColumn){
        return M[ithRow][jthColumn];
    }
    Matrix(int rows = 0, int columns = 0){
        write(rows, columns);
    }
    Matrix(int rows, int columns, T** array){
        write(rows, columns);
        take_array(array);
    }
    void take_array(T** array){
        for(int i = 0; i < m; i++)
        for(int j = 0; j < n; j++)
        M[i][j] = array[i][j];
    }
    bool is_empty(){
        return !m;
    }
    bool is_square(){
        return (!is_empty() && m == n);
    }
    bool is_singlet(){
        return (is_square() && m == 1);
    }
    int row_lagging_zeros(int i){
        if(!is_empty()){
            int count = 0;
            for(int j = 0; j < n; j++){
                if(M[i][j])return count;
                count++;
            }
            return -2;
        }
        else return -1;
    }
    T row_leading_entry(int i){
        if(!is_empty()){
            for(int j = 0; j < n; j++){
                if(M[i][j])return M[i][j];
            }
            return 0;
        }
        return numeric_limits<T>::min();
    }
    bool is_zero_row(int i){
        if(!is_empty()){
            bool yes = true;
            for(int j = 0; j < n; j++){
                if(M[i][j])return false;
            }
            return yes;
        }
        else return false;
    }
    void swap_rows(int i0, int i1){
        if(is_empty())return;
        for(int j = 0; j < n; j++){
            swap(M[i0][j], M[i1][j]);
        }
    }
    void pivot_rows(int i0, int i1, T multiple){
        if(is_empty())return;
        for(int j = 0; j < n; j++){
            M[i0][j] += M[i1][j] * multiple;
        }
    }
    void rescale_rows(int i0, T scale){
        if(is_empty())return;
        for(int j = 0; j < n; j++){
            M[i0][j] *= scale;
        }
    }
    void swap_columns(int j0, int j1){
        if(is_empty())return;
        for(int i = 0; i < m; i++){
            swap(M[i][j0], M[i][j1]);
        }
    }
    void pivot_columns(int j0, int j1, T multiple){
        if(is_empty())return;
        for(int i = 0; i < m; i++){
            M[i][j0] += M[i][j1] * multiple;
        }
    }
    void rescale_columns(int j0, T scale){
        if(is_empty())return;
        for(int i = 0; i < m; i++){
            M[i][j0] *= scale;
        }
    }
    static Matrix Multiply(Matrix mat1, Matrix mat2){
        if(mat1.is_empty()||mat2.is_empty())return Matrix();
        if(mat1.n == mat2.m){
            Matrix product(mat1.m, mat2.n);
            for(int i = 0; i < mat1.m ; i++){
                for(int j = 0; j < mat2.n; j++){
                    T dot_product = 0;
                    for(int k = 0; k < mat2.m; k++){
                        dot_product += mat1.M[i][k] * mat2.M[k][j];
                    }
                    product.M[i][j] = dot_product;
                }
            }
            return product;
        }
        else return Matrix();
    }
    Matrix Inverse(){
        if(!is_empty() && is_square()){
            T det = determinant();
            if(det){
                return adjacent() * (1 / det);
            }
            else return Matrix();
        }
        else return Matrix();
    }
    bool is_row_reduced_echelon(){
        if(!is_empty()){
            if(!is_row_echelon())return false;
            for(int i = 0; i < m; i++){
                for(int j = 0; j < n; j++){
                    if(M[i][j] == 1){
                        for(int k = 0; k < m; k++){
                            if(i == k)continue;
                            else if(M[k][j])return false;
                        }
                    }
                }
            }
            return true;
        }
        else return false;
    }
    bool is_zero_column(int j){
        if(!is_empty()){
            bool yes = true;
            for(int i = 0; i < m; i++){
                if(M[i][j])return false;
            }
            return yes;
        }
        else return false;
    }
    bool is_row_echelon(){
        if(!is_empty()){
            for(int i = 0; i < m; i++){
                if(row_leading_entry(i) != 1)
                {
                    if(i != m - 1)return false;
                    else if(!is_zero_row(i))return false;
                }
                else if(i && row_lagging_zeros(i) <= row_lagging_zeros(i - 1)){return false;}
            }
            return true;
        }
        else
        {
            return false;
        } 
    }
    Matrix adjacent(){
        return cofactor().transpose();
    }
    Matrix transpose(){
        if(!is_empty()){
            Matrix trans(n, m);
            for(int i = 0; i < n; i++){
                for(int j = 0; j < m; j++)
                trans.M[i][j] = M[j][i];
            }
            return trans;
        }
        else return Matrix();
    }
    T cofactor(int i, int j){
        return (is_square() ? ((((i + j) % 2) ? -1 : 1) * minor(i, j).determinant()) : numeric_limits<T>::min());
    }
    Matrix cofactor(){
        if(is_square()){
            Matrix mat(m, m);
            for(int i = 0; i < m; i++){
                for(int j = 0; j < m; j++){
                    mat.M[i][j] = (((i + j) % 2) ? -1 : 1) * minor(i, j).determinant();
                }
            }
            return mat;
        }
        else return Matrix();
    }
    Matrix minor(int k, int l){
        if(is_square()){
            if(is_singlet())return *this;
            T** Matn_1 = new T*[m - 1];
            for(int i = 0; i < m - 1; i++)
            Matn_1[i] = new T[m - 1];
            int ic = 0;
            for(int i = 0; i < m; i++){
                int jc = 0;
                for(int j = 0; j < m; j++){
                    if(i != k && j != l){
                        Matn_1[ic][jc++] = M[i][j];
                    }
                }
                if(i != k)ic++;
            }
            return Matrix(m - 1, m - 1, Matn_1);
        }
        else return Matrix();
    }
    T determinant(){
        if(!is_square())return numeric_limits<T>::min();
        else {
            T det = 0;
            if(m == 1) det = M[0][0];
            else
            for(int j = 0; j < m; j++)
            {
                det += (j % 2 ? -1 : 1) * M[0][j] * minor(0, j).determinant();
            }
            return det;
        }
    }
    void write(int rows = 0, int columns = 0){
        m = rows;
        n = columns ? (rows ? columns : rows) : rows;
        M = (rows ? new T*[rows] : nullptr);
        for(int i = 0; i < rows; i++){
            M[i] = new T[columns];
        }
    }
    void rewrite(int rows = 0, int columns = 0){
        for(int i = 0; i < m; i++){
            delete M[i];
        }
        delete[] M;
        write(rows, columns);
    }
    T** getMatrix() {
        return M;
    }
    ushort getColumns() const{
        return n;
    }
    ushort getRows() const{
        return m;
    }
    template <class Z>
    friend Matrix<Z> operator*(Matrix<Z>, Matrix<Z>);
    template <class Z>
    friend Matrix<Z> operator*(Matrix<Z>, Z);
    Matrix operator=(Matrix mat){
        rewrite(mat.m, mat.n);
        make_copy(mat);
        return *this;
    }
    friend istream& operator>>(istream&, Matrix<float>&);
    template <class Z>
    friend ostream& operator<<(ostream&, Matrix<Z>);
};
template <class Z>
Matrix<Z> operator*(Matrix<Z> mat1, Matrix<Z> mat2){
    return Matrix<Z>::Multiply(mat1, mat2);
}
template <class Z>
Matrix<Z> operator*(Matrix<Z> mat, Z scale){
    Matrix<Z> mat2 = mat;
    for(int i = 0; i < mat2.m ; i++){
        mat2.rescale_rows(i, scale);
    }
    return mat2;
}
istream& operator>>(istream& in, Matrix<float>& mat){
    if(!mat.getColumns() || !mat.getRows()){
        ushort x, y;
        cout << " NULL MATRIX!!!\a\n";
        cout << " Enter size (rows x columns): ";
        in >> x >> y;
        mat.rewrite(x, y);
    }
    cout << " Enter " << mat.getRows() <<" x "<< mat.getColumns() << " : \n" ;
    for(int i = 0; i < mat.getRows(); i++){
        cout << ' ' << i + 1 <<"th row: ";
        for(int j = 0; j < mat.getColumns(); j++){
            in >> mat.giveEntry(i, j);
        }
    }
    return in;
}
template <class Y>
ostream& operator<<(ostream& out, Matrix<Y> mat){
    if(!mat.getColumns() || !mat.getRows()){
        out << "EMPTY MATRIX";
        return out;
    }
    for(int i = 0; i < mat.getRows(); i++){
        for(int j = 0; j < mat.getColumns(); j++){
            out << mat.getMatrix()[i][j] << ' ';
        }
        out << endl ;
    }
    return out;
}
ifstream& operator>>(ifstream& ifs, Matrix<float>& mat){
    if(!mat.getColumns() || !mat.getRows()){
        ushort x, y;
        cout << " NULL MATRIX!!!\a\n";
        cout << " Enter size (rows x columns): ";
        cin >> x >> y;
        mat.rewrite(x, y);
    }
    for(int i = 0; i < mat.getRows(); i++)
    {
        for(int j = 0; j < mat.getColumns(); j++)
        {
            ifs >> mat.getMatrix()[i][j];
        }
    }
    return ifs;
}