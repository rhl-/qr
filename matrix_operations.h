#ifndef T10_MATRIX_OPERATIONS
#define T10_MATRIX_OPERATIONS
#include <iostream>
#include <math.h>
#include <stdio.h>

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
	if (s == 0.0) {
		b = 0.0;
	}
	else {
		double mu = sqrt(x[0]*x[0] + s);
		if (x[0] <= 0.0) {
			v[0] = x[0] - mu;
		}
		else {
			printf("s=%f\n",s);
			printf("x0=%f\n",x[0]);
			printf("mu=%f\n",mu);
			v[0] = -(double)s/(double)(x[0]+mu);
			printf ("v[0]=%f\n", v[0]);
		}
		b = 2.0*v[0]*(double)v[0]/(double)(s + v[0]*v[0]);
		printf ("v[0]=%f\n",v[0]);
		for (std::size_t i = 0; i < v.size(); ++i) {
			printf ("v[i]=%f\n",v[i]);
			v[i] = (double)v[i]/(double)v[0];
		}
	}
}
 //end namespace T10 namespace
}
#endif //ifndef T10_MATRIX_OPERATIONS
