#ifndef T10_MATRIX_OPERATIONS
#define T10_MATRIX_OPERATIONS
#include <iostream>
namespace t10 {

template<class Matrix, class Vector> 
void mat_vec_mult (const Matrix & A, const Vector & b, Vector & c) {
	for (std::size_t i = 0; i < c.size(); ++i){
		c[i] = 0.0;
		for( std::size_t j = 0; j < b.size(); ++j){
			c[i] += A[i][j]*b[j];
		}
	}	
}

} //end namespace T10 namespace
#endif //ifndef T10_MATRIX_OPERATIONS
