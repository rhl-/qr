#ifndef T10_MATRIX_OPERATIONS
#define T10_MATRIX_OPERATIONS
#include <iostream>
#include <math.h>

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

template<class Vector>
double house(Vector & v) {
	const std::size_t n = v.size();
	double beta = 0.0;
	double sigma = 0.0;
	for (std::size_t i = 1; i < v.size(); ++i) { sigma += v[i]*v[i]; };
	const double x = v[0];
	if (sigma != 0){
		const double mu = sqrt(x*x + sigma);
		(x <= 0)? v[0] -= mu: v[0] = -sigma/(x+mu);
		const double y = v[0]*v[0];
		const beta = (2.0*y)/(sigma + y);
		for (std::size_t i = 1; i < v.size(); ++i) { v[i] = v[i]/v[0]; }
		v[0] = 1;
	}
	return beta;
}
 //end namespace T10 namespace
}
#endif //ifndef T10_MATRIX_OPERATIONS
