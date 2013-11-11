#include <iostream>
#include <vector>
#include <iomanip> 
#include "matrix_operations.h"

typedef std::vector<double> Vector;
int main () {
std::vector< Vector > A(3,Vector(2,1));

std::cout << std::fixed << std::setprecision( 5 );

Vector b;
b.push_back(3.0);
b.push_back(2.0);
b.push_back(1.0);

Vector c(2,0);
t10::mat_vec_mult(A, b, c);
for (int i = 0; i < A.size(); i++){
	for (int j = 0; j < A[0].size(); j++){
		std::cout << A[i][j] << " ";
	}
	std::cout << std::endl;
}

Vector v(3,0.0);
double beta;
t10::house(b,v,beta);
std::cout << "beta = " << beta << std::endl;
for (typename Vector::const_iterator i  = v.begin(); i != v.end(); ++i) {
	std::cout << *i << std::endl;
} 
return 0;
}

