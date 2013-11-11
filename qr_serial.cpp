
#include <iostream>
#include "qr_algorithm.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include<algorithm>

namespace ublas = boost::numeric::ublas;

typedef ublas::vector< double > Vector;
typedef ublas::matrix< double> Matrix;
typedef ublas::zero_matrix< double> Zero_matrix;
typedef ublas::zero_vector< double> Zero_vector;
int main( int argc, char * argv[]){
	const double test_case[3][3] = 
       {{5,3,3},
	{8,7,4},
	{4,3,8}};
	 
	const double correct_hess[3][3] = 
  	{{ 5.000000000000000,   -4.024922359499620,    1.341640786499874},
   	 {-8.944271909999159,   10.000000000000000,   -3.000000000000000},
   	 { 0.000000000000000,   -2.000000000000000,    5.000000000000000}};

	const double vector_test[3] = {3,2,1};
   	const double correct_house[3] = {1.00000000000000,-2.69666295470958,-1.34833147735479};
	const double correct_beta =  0.198216274262727;

	//define matrix type (use boost)
	//read input with boost program options
	//boost supports specializations for banded, triangular, etc, see:
	//http://www.boost.org/doc/libs/1_54_0/libs/numeric/ublas/doc/index.htm
	const std::size_t n = 3;
	Matrix M(n,n);
	Matrix A(n,n);
	Vector V(n);
	Vector CV(n);
	for (unsigned int i = 0; i < M.size1(); ++i){
		V(i) = vector_test[i];
		CV(i) = correct_house[i];
		for(unsigned int j = 0; j < M.size2(); ++j){
			M(i,j) = test_case[i][j];
			A(i,j) = correct_hess[i][j];
		}
	}
	//Householder Vector Test
	std::cout << "Input: " << V << std::endl;
	t10::compute_householder_vector(V);
	double beta = V(0); V(0) = 1;
	std::cout << "beta: " << beta << " correct: " << correct_beta << std::endl;
	std::cout << "V: " << V << " correct V: " << CV << std::endl;
	std::cout << " ------------------------------------------------- " << std::endl;	
	//Hessenberg Reduction Test
	std::cout << "Input Matrix: " << M << std::endl;
	t10::qr(M);
	std::cout << "M = " << M << std::endl;
	std::cout << "A = " << A << std::endl;
}