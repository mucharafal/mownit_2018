//#include <boost/multiprecision/cpp_dec_float.hpp> //not so easy to use
#include <iostream>
#include <iomanip>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

template <class T>
class denseMatrix{
  public:
    T** matrix;
    int size;
    denseMatrix(int size){
        this->size = size;
        this->matrix = (T**)malloc(size * sizeof(T*));
        for(int i=0;i< size;i++){
            this->matrix[i] = (T*)malloc(size * sizeof(T));
        }
    }
    ~denseMatrix(){
        for(int i = 0;i < size;i++){
            free(matrix[i]);
        }
        free(matrix);
    }
    void print(){
        for(int i =0;i < size;i++){
            for(int j = 0;j < size;j++)
                cout << matrix[i][j] << " ";
            cout << endl;
        }
    }

    void generate(int case_number){
        switch(case_number){
            case 1:
            for(int i = 0;i < size;i++){
                matrix[0][i] = 1;
                for(int j = 1;j < size;j++){
                    matrix[j][i] = ((T)1) / (i + j + 1);
                }
            }
            break;
            case 2:
            for(int i = 0;i < size;i++){
                matrix[i][i] = 2;
                for(int j = i+1;j < size;j++){
                    matrix[j][i] = matrix[i][j] = ((T)2) * (i+1) / (j+1);
                }
            }
            break;
            case 3:
            int k = 7, m = 3;
            matrix[0][0] = k;
            matrix[0][1] = ((T)1) / (1+m);
            for(int i=1;i<(size-1);i++){
                matrix[i][i] = k;
                matrix[i][i+1] = ((T)1) / (i+1+m);
                matrix[i][i-1] = (T(k)) / (i+m+2);
            }
            matrix[size-1][size-1] = k;
            matrix[size-1][size-2] = (T(k)) / ((size-1)+m+2);
        }
    }
};

int main(){
    int size;
    int triesNumber, type;
    cin >> triesNumber;
    for(int i = 0;i < triesNumber;i++ ){
        cin >> size >> type;
        mat<double> A(size, size);
        denseMatrix<double>* matrix = new denseMatrix<double>(size);
        matrix->generate(type);
        for(int i = 0;i < size;i++){
            for(int j = 0;j < size;j++)
                A(i, j) = matrix->matrix[i][j];
        }
        cout << cond(A) << endl;
    }
}
