#include <boost/cstdfloat.hpp> // For float_64_t, float128_t. Must be first include!
#include <boost/config.hpp>
#include <boost/multiprecision/float128.hpp>
//#include <boost/multiprecision/cpp_dec_float.hpp> //not so easy to use
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace boost::multiprecision;

template <class T>
T* substractVectors(T* vecA, T* vecB, int size){
    T* result = (T*)malloc(size * sizeof(T));
    for(int i = 0;i < size;i++){
        result[i] = vecA[i] - vecB[i];
    }
    return result;
}

template <class T>
T multipleVectors(T* vecA, T* vecB, int size){
        T res = 0;
        for(int i = 0;i < size;i++) res += vecA[i] * vecB[i];
        return res;
}

template <class T>
T euklidesNorm(T* vec, int size){
    return (T)sqrt(multipleVectors(vec, vec, size));
}

template <class T>
T maksimumNorm(T* vec, int size){
    T maks = 0;
    for(int i = 0;i < size;i++){
        maks = max(maks, vec[i]);
    }
    return maks;
}

template <class T>
T* generateOnesVector(int size){
    T* vectorB = (T*)malloc(sizeof(T) * size);
    for(int i = 0;i < size;i++){
        if(i % 2 == 1){
            vectorB[i] = 1;
        } else {
            vectorB[i] = -1;
        }
    }
    return vectorB;
}

template <class T>
class matrixClass{
    public:
    virtual T* multipleVector(T* vectorX) = 0;
    virtual T* generate(int, int) = 0;
    virtual T* solve(T*) = 0;
};

template <class T>
class denseMatrix: public matrixClass<T>{
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

    void swapTwoRows(int first, int second){
        swap(matrix[first], matrix[second]);
    }

    void minus(int from, int a, T aTimes){
        T* r1 = matrix[from];
        T* r2 = matrix[a];

        for(int i = 0;i < size;i++){
            r1[i] -= r2[i] * aTimes;
        }
    }

    T* multipleVector(T* vector){
        T* result = (T*)malloc(sizeof(T) * this->size);
        for(int i = 0;i < this->size;i++)    result[i] = 0;
        for(int row=0;row < this->size;row++){
            for(int column=0;column < size;column++){
                result[row] += (this->matrix[row][column] * vector[column]);
            }
        }
        return result;
    }

    void prepareGaussianMatrix(T* vectorB){
        for(int column = 0;column < size;column++){

            int maximumIndex = findMaximumInSubmatrix(column);
            swapTwoRows(column, maximumIndex);
            swap(vectorB[column], vectorB[maximumIndex]);

            for(int row = column + 1;row < this->size;row++){
                T diff = this->matrix[row][column] / this->matrix[column][column];
                minus(row, column, diff);
                vectorB[row] -= vectorB[column] * diff;
            }
        }
    }

    int findMaximumInSubmatrix(int startIndex){
        T maximum = 1;
        int maximumIndex = startIndex;
        for(int i = startIndex;i < size;i++){
            for(int j = startIndex;j < size;j++){
                T value = abs(matrix[i][startIndex] / matrix[i][j]);
                if(maximum < value){
                    maximum = value;
                    maximumIndex = i;
                }
            }
        }
        return maximumIndex;
    }

    T* solve(T* vectorB){
        prepareGaussianMatrix(vectorB);
        T* x = subsolve(vectorB);
        return x;
    }

    T* subsolve(T* vectorB){
        T* x = (T*)malloc(sizeof(T) * size);
        for(int i = 0;i < size;i++)    x[i] = 0;
        for(int row = size-1;row >= 0;row --){
        //T[row][row] * x + a = b
            T a = multipleVectors(matrix[row], x, size);
            x[row] = (vectorB[row] - a) / matrix[row][row];
        }
        return x;
    }

    void print(){
        for(int i =0;i < size;i++){
            for(int j = 0;j < size;j++)
                cout << matrix[i][j] << " ";
            cout << endl;
        }
    }
    template <class R>
    void print(R* vecB){
        for(int i =0;i < size;i++){
            for(int j = 0;j < size;j++)
                cout << matrix[i][j] << " ";
            cout << vecB[i] << endl;
        }
    }

    T* generate(int case_number, int junk){
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
        return generateOnesVector<T>(size);
    }
};

template <class T>
class diagonalMatrix: public matrixClass<T> {
  private:
    T* upperDiagonal;
    T* mediumDiagonal;
    T* lowerDiagonal;
    int size;
  public:
    diagonalMatrix(int size){
        upperDiagonal = (T*)malloc(sizeof(T) * (size-1));
        lowerDiagonal = (T*)malloc(sizeof(T) * (size-1));
        mediumDiagonal = (T*)malloc(sizeof(T) * (size));
        this->size = size;
    }
    ~diagonalMatrix(){
        free(upperDiagonal);
        free(lowerDiagonal);
        free(mediumDiagonal);
    }

    T* generate(int o, int w){
        int k = 7;
        int m = 3;
        int ssize = size-1;
        for(int i=0;i<ssize;i++){
            mediumDiagonal[i] = k;
            upperDiagonal[i] = ((T)1) / (i+1+m);
            lowerDiagonal[i] = (T(k)) / (i+m+2);
        }
        mediumDiagonal[size-1] = k;

        return generateOnesVector<T>(size);
    }

    void prepareGaussianMatrix(T* vectorB){
        for(int i=0;i<(size-1);i++){
            T multiplier = lowerDiagonal[i] / mediumDiagonal[i];
            lowerDiagonal[i] = 0;
            mediumDiagonal[i+1] -= upperDiagonal[i] * multiplier;
            vectorB[i+1] -= vectorB[i] * multiplier;
        }
    }

    T* subsolve(T* vectorB){
        T* x = (T*)malloc(size * sizeof(T));
        x[size-1] = vectorB[size-1] / mediumDiagonal[size-1];
        for(int i=size-2;i>=0;i--){
            x[i] = (vectorB[i] - upperDiagonal[i] * x[i+1]) / mediumDiagonal[i];
        }
        return x;
    }

    T* solve(T* vectorB){
        prepareGaussianMatrix(vectorB);
        return subsolve(vectorB);
    }

    T* multipleVector(T* vectorX){
        T* result = (T*)malloc(size * sizeof(T));
        result[0] = mediumDiagonal[0] * vectorX[0] + upperDiagonal[0] * vectorX[1];
        for(int i = 1;i < (size-1);i++){
            result[i] = lowerDiagonal[i-1] * vectorX[i-1] +
                mediumDiagonal[i] * vectorX[i] +
                upperDiagonal[i] * vectorX[i+1];
        }
        result[size-1] = lowerDiagonal[size-2] * vectorX[size-2] +
            mediumDiagonal[size-1] * vectorX[size-1];
        return result;
    }
};

template <class T>
void makeTest(int matrixType, int method, int size){
    matrixClass<T>* matrix;
    if(matrixType == 1)
        matrix = new diagonalMatrix<T>(size);
    else
        matrix = new denseMatrix<T>(size);
    T* vecX = matrix->generate(method, 3);

    T* vecB = matrix->multipleVector(vecX);
    clock_t start_time = clock();
    T* newX = matrix->solve(vecB);
    clock_t end_time = clock();

    T* diff = substractVectors(vecX, newX, size);



    cout << " ; "<< euklidesNorm(diff, size) << " ; " << maksimumNorm(diff, size) << " ; " << (end_time - start_time) / (double)(CLOCKS_PER_SEC / 1000) << endl;

    free(vecX);
    free(newX);
    free(diff);
    if(matrixType == 1)
        delete (diagonalMatrix<T>*)matrix;
    else
        delete (denseMatrix<T>*)matrix;
}

int main()
{
    int caseNumber;
    cin >> caseNumber;
    int matrixType;
    int generateType;
    int precision;
    cin >> matrixType >> generateType >> precision;
    for(int i = 0;i < caseNumber;i++){
        int size;
        cin >> size;

        cout << size << "  ";
        switch(precision){
            case 1:
            makeTest<float>(matrixType, generateType, size);
            break;
            case 2:
            makeTest<double>(matrixType, generateType, size);
            break;
			case 3:
			makeTest<float128>(matrixType, generateType, size);
			break;
            case 4:
            makeTest<long double>(matrixType, generateType, size);
            break;
        }
    }
    return 0;
}
