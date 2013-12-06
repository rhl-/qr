
#include <iostream>
#include "qr_serial.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <algorithm>
#include <boost/numeric/ublas/banded.hpp> //for banded_adaptor

namespace ublas = boost::numeric::ublas;
namespace serial = t10::serial;

typedef ublas::vector< double > Vector;
typedef ublas::matrix< double> Matrix;
typedef ublas::zero_matrix< double> Zero_matrix;
typedef ublas::zero_vector< double> Zero_vector;
typedef typename ublas::diagonal_adaptor< Matrix> Diagonal_adapter;
int main( int argc, char * argv[]){
	const double test_case[3][3] = 
       {{5,3,3},
	{8,7,4},
	{4,3,8}};
	Matrix M0(0,0);
	Matrix M1(1,1);
	M1(0,0) = 42;
	Matrix M2(2,2);
	M2(0,0) = 42;  
	M2(0,1) = 24; 
	M2(1,0) = 69; 
	M2(1,1) = 12;
	serial::qr(M1);
	serial::qr(M2);
	std::cout << t10::print_matrix(M1) << std::endl << t10::print_matrix(M2) << std::endl;
	const double correct_hess[3][3] = 
  	{{ 5.000000000000000,   -4.024922359499620,    1.341640786499874},
   	 {-8.944271909999159,   10.000000000000000,   -3.000000000000000},
   	 { 0.000000000000000,   -2.000000000000000,    5.000000000000000}};

	const double vector_test[3] = {3,2,1};
   	const double correct_house[3] = {1.00000000000000,-2.69666295470958,-1.34833147735479};
	const double correct_beta =  0.198216274262727;
	const double correct_eigs[3] = {14.6235, 1.0000, 4.3765};
	//define matrix type (use boost)
	//read input with boost program options
	//boost supports specializations for banded, triangular, etc, see:
	//http://www.boost.org/doc/libs/1_54_0/libs/numeric/ublas/doc/index.htm
	const std::size_t n = 3;
	Matrix M(n,n);
	Matrix HM(n,n);
	Matrix A(n,n);
	Vector V(n);
	Vector CV(n);
	Vector E(n);
	for (unsigned int i = 0; i < M.size1(); ++i){
		V(i) = vector_test[i];
		CV(i) = correct_house[i];
		E(i) = correct_eigs[i];
		for(unsigned int j = 0; j < M.size2(); ++j){
			M(i,j) = test_case[i][j];
			A(i,j) = correct_hess[i][j];
		}
	}
	HM = M;
	//Householder Vector Test
	std::cout << "Input: " << V << std::endl;
	t10::serial::compute_householder_vector(V);
	double beta = V(0); V(0) = 1;
	std::cout << "beta: " << beta << " correct: " << correct_beta << std::endl;
	std::cout << "V: " << V << " correct V: " << CV << std::endl;
	std::cout << " ------------------------------------------------- " << std::endl;	
	//Hessenberg Reduction Test
	std::cout << "Input Matrix: " << t10::print_matrix(HM)<< std::endl;
	t10::serial::hessenberg(HM);
	std::cout << "M = " << t10::print_matrix(HM) << std::endl;
	std::cout << "A = " << t10::print_matrix(A) << std::endl;
	t10::serial::qr(M);
	Diagonal_adapter D(M);
	std::cout << "D = " << t10::print_matrix(D) << std::endl;
	std::cout << "Correct Eigs: " << E << std::endl; 	
}
