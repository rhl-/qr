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
	double s = 0.0;
	v[0] = 1.0;
	for (std::size_t i = 1; i < x.size(); ++i) {
		s += x[i]*x[i];
		v[i] = x[i];
	}
	std::cout << s << std::endl;
	std::cout << v[1] << std::endl;
	if (s == 0.0) { b = 0.0; }
	else {
		double mu = sqrt(x[0]*x[0] + s);
		if (x[0] <= 0.0) {
			v[0] = x[0] - mu;
		}
		else {
			std::cout << "s " << s << std::endl;
			std::cout << "x0 " << x[0] << std::endl;
			std::cout << "mu " << mu << std::endl;
			v[0] = -s/(x[0]+mu);
			std::cout << "v[0] " << v[0] << std::endl;
		}
		b = (2.0*v[0]*v[0])/(s + v[0]*v[0]);
		std::cout << "v[0] " << v[0] << std::endl;
		for (std::size_t i = 0; i < v.size(); ++i) {
			v[i] = (double)v[i]/(double)v[0];
			std::cout << "v[i] " << v[i] << std::endl;
		}
	}
}
 //end namespace T10 namespace
}
#endif //ifndef T10_MATRIX_OPERATIONS
