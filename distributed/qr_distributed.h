#ifndef QR_DISTRIBUTED_H
#define QR_DISTRIBUTED_H
#define QR_ITERATION_OUTPUT

//BOOST MPI
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

//BOOST UBLAS
#include <boost/numeric/ublas/banded.hpp> //for banded_adaptor
#include <boost/numeric/ublas/matrix.hpp> //for slices
#include <boost/numeric/ublas/matrix_proxy.hpp> //for slices
#include <boost/numeric/ublas/vector_proxy.hpp> //for slices

//PROJECT
#include "io.h"
#include "qr_serial.h"
#include "hessenberg_distributed.h"

namespace ublas = boost::numeric::ublas;
namespace t10 {
	namespace parallel {
	template< typename Matrix>
	typename Matrix::value_type compute_shift( const Matrix & H){
		typedef typename Matrix::value_type Value;
		const std::size_t n = H.size1();
		const Value a = H(n-2,n-2)+H(n-1,n-1);
		const Value b = std::sqrt(4*H(n-1,n-2)*H(n-2,n-1) + 
				std::pow((H(n-2,n-2)-H(n-1,n-1)),2));
		const Value e1 = (a+b)/2, e2 = (a-b)/2;
		const Value l = std::fabs(e1- H(n-1,n-1));
		const Value r = std::fabs(e2 - H(n-1,n-1));
		return (l<r)? e1:e2;
	}
	template< typename Matrix_data>
	void shift_matrix( Matrix_data & data, int sign=-1){
		typedef typename Matrix_data::Matrix Matrix;
		typedef typename Matrix::value_type Value;
		typedef typename ublas::diagonal_adaptor< Matrix> 
						  Diagonal_adaptor;
		typedef typename ublas::scalar_matrix< Value> Scalar;
		typedef ublas::diagonal_adaptor< Scalar> Diagonal_scalar;
		Matrix & H = data.M;
		const std::size_t n = H.size1();
		const std::size_t p = data.row_length;
		Value shift = compute_shift( H);
		mpi::broadcast( data.diag_comm, shift, p-1);
		Scalar mu_s(n,n,shift);
		const Diagonal_scalar mu( mu_s );
		Diagonal_adaptor D(H);
		D += sign*mu;
	}

	template< typename Matrix_data>
	void qr_iteration( Matrix_data & data, double tol=1e-16){
		typedef typename Matrix_data::Matrix Matrix;	
		typedef typename Matrix::value_type Value;
		Matrix & H = data.M;
		#ifdef QR_ITERATION_OUTPUT
		std::cout << "H = " << t10::print_matrix( H) << std::endl;
		#endif
		const std::size_t n = H.size1();
		std::vector< typename Matrix::value_type> givens(n-1, 0.0); 
		bool done = false;
		do{
			if (data.diag()){ shift_matrix( data); }
			//Step 0: Ghost boundary rows and columns
			//Step 1: (Diagonals) 
				//a) Wait for signal to start,
				//b) apply_g_left, to build givens matrices
				//c)  bcast data: 
					//1) across row_comm[ row_index]
					//2) up col_comm[ row_index]
				//d) apply_g_right
			//Step 1: (Off Diagonals)
				//NOTE: each of (a) and (b) bcasts are 
				//a) bcast for data from row diagonal processor
				//   if( non_empty){ apply_g_left }
				//   else { break; }
			
				//b) bcast for data from diagonal col processor
				//   apply_g_right
			//Step 2: (Diagonals)
			
			for (std::size_t i = 0; i < n-1; ++i){ 
				givens[i]=apply_givens_left(H,i,i+1);
			}
			
			#ifdef DEBUG_QR_ITERATION
			std::cout << "Left Apply: H = " 
				  << t10::print_matrix( H) << std::endl;
			#endif //DEBUG_QR_ITERATION
			
			for (std::size_t i = 0; i < n-1; ++i){ 
				apply_givens_right(H,givens[i],i,i+1); 
			}
			
			#ifdef DEBUG_QR_ITERATION
			std::cout << "Right Apply H = " 
				  << t10::print_matrix( H) << std::endl;
			#endif //DEBUG_QR_ITERATION
			
			//Unshift
			if( data.diag()){ shift_matrix( data, 1); }
			
			#ifdef QR_ITERATION_OUTPUT
			//std::cout.precision( 7);
			//std::cout.setf(std::ios::fixed,std:: ios::floatfield);
			//std::cout << "---- Shift: " << std::setw( 10) 
			//	  << shift  << std::endl;
			//std::cout << "---- Error: " << std::setw( 10) 
			//	  << std::fabs(H(n-1,n-2)) << std::endl;
			#endif //DEBUG_QR_ITERATION
			if( data.diag()){
				//std::size_t idx = in_charge(data, i); 
				//mpi::broadcast( data.diag_comm, done, idx); 
				if (done){
					//bcast empty givens package across row
				}
			}
		} while( !done); 
	}
	
	template< typename Matrix_data>
	void qr( Matrix_data & data){
			hessenberg( data);
			//for (std::size_t i = data.n; i > 1; --i){
			//	Range r( begin, end);
			//	Matrix_range R(M, r, r);
			//qr_iteration( data);
			//}
	}

} //end namespace parallel
} //end namespace t10
#endif //QR_DISTRIBUTED_H
