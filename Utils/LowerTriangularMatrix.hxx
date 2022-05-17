/*
	This module contains efficient structures to represent.

	Bobak Pezeshki, May 2022.
*/

#ifndef LOWER_TRIANGULAR_MATRIX_HXX_INCLUDED
#define LOWER_TRIANGULAR_MATRIX_HXX_INCLUDED

#include <inttypes.h>

#ifdef LINUX
typedef int64_t INT64 ;
#else 
typedef __int64 INT64 ;
#endif

// lower triangular matrix without diagonal
template<typename T>
class LowerTriangularMatrix {
    protected:
        T* arr;
        uint32_t d;

    public:
        LowerTriangularMatrix(uint32_t dim_) {
            arr = new T[(dim_ - 1)*dim_/2];
            d = dim_;
        }

    public:
        T &operator()(uint32_t i, uint32_t j) {
            // idx(i,j) = (d-1)*d/2 - (d-i-1)*(d-i)/2 + j - i - 1 = (2*d*i - i*i - i)/2 + j - i - 1

            // ex. d=5
            // - 0 0 0 0
            // - - 0 0 0
            // - - - 0 0
            // - - - - 0
            // - - - - -

            // [0][1] = 0
            // [0][2] = 1
            // [0][3] = 2
            // [0][4] = 3
            // [1][2] = 4  (4)*5/2 - (3)*(4)/2 + 2 - (1) - 1 = 10 - 6 + 2 - 1 - 1 = 4
            // [1][3] = 5
            // [1][4] = 6  (4)*5/2 - (3)*(4)/2 + 4 - (1) - 1 = 10 - 6 + 4 - 1 - 1 = 6
            // [2][3] = 7
            // [2][4] = 8  (5-1)*5/2 - (5-2-1)*(5-2)/2 + 4 - 2 - 1 = 10 - 3 + 4 - 2 - 1 = 8
            // [3][4] = 9

            return arr[(2*d*i - i*i - i)/2 + j - i - 1];
        };

        uint32_t len() const{
            return (d - 1)*d/2;
        }

        T* arr_ptr() {
            return arr;
        }

        T val(uint32_t i, uint32_t j) const {
            return arr[(2*d*i - i*i - i)/2 + j - i - 1];
        };

        ~LowerTriangularMatrix(){
            delete[] arr;
        }
};


#endif // LOWER_TRIANGULAR_MATRIX_HXX_INCLUDED
