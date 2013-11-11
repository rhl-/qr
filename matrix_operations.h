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
void house(const Vector & x, Vector & v, double & b) {
	double y;
	double s = 0.0;
	b = 0.0;
	for (std::size_t i = 1; i < x.size(); ++i) {
		s += x[i]*x[i];
		v[i] = x[i];
	}
	v[0] = 1.0;
	if (s > 0) {
		const double mu = sqrt(x[0]*x[0] + s);
		(x[0] <= 0.0) ? v[0] = x[0] - mu : v[0] = -s/(x[0]+mu);
		const double v0 = v[0];
		y = v0*v0;
		b = 2.0*y/(s + y);
		for (std::size_t i = 0; i < v.size(); ++i) {
			v[i] = v[i]/v0;
		}
	}
}
 //end namespace T10 namespace
}
#endif //ifndef T10_MATRIX_OPERATIONS
