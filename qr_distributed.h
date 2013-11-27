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

namespace ublas = boost::numeric::ublas;
namespace t10 {
	namespace parallel {
	//GVL Section 5.1.9
	template< typename Value>
	void compute_givens(const Value & a, const Value & b, 
			    Value & c, Value & s){
			if (b == 0) { c = 1; s = 0; return; }
			if( std::fabs(a) < std::fabs(b)){
				const Value tau = -a/b;
				s = 1/ std::sqrt(1 + std::pow(tau, 2));
				c = tau * s;
				return;
			}
			const Value tau = -b/a;
			c = 1/ std::sqrt(1 + std::pow(tau, 2));
			s = tau * c;
	}

	template< typename T>
	int sign( const T & x){ return (x < 0) ? -1 : 1; }
	
	template< typename Value>
	Value encode_givens( const Value & c, const Value & s){
			if (c == 0) { return 1; }
			if (std::fabs(s) < std::fabs(c)) { return sign(c)*s/2; }
			return 2*sign(s)/c;
	}
	template< typename Value>
	void decode_givens( const Value & rho, Value & c, Value & s ){
		if (rho == 1) { c = 0; s = 1; return; }
		if (std::fabs(rho) < 1) { 
			s = 2*rho; 
			c = std::sqrt(1 -std::pow(s,2)); 
			return;
		}
		c = 2/rho; s = std::sqrt(1-std::pow(c,2));
	} 
	
	template< typename Matrix>
	typename Matrix::value_type apply_givens_left( Matrix & M, 
						       const std::size_t i, 
						       const std::size_t k){
		typedef typename Matrix::value_type Value;
		Value c=0.0, s=0.0;
		compute_givens(M(i,i),M(k,i),c,s);
		for(std::size_t j = i; j < M.size2(); ++j){
			const double ti = M(i,j);
			const double tj = M(k,j);
			M(i,j) = c*ti - s*tj;
			M(k,j) = s*ti + c*tj;
		}
		M(k,i) = 0;
		return encode_givens(c,s);
	}
	template< typename Matrix>
	void apply_givens_right( Matrix & M, 
				 const typename Matrix::value_type rho, 
				 const std::size_t i, const std::size_t k){
		typename Matrix::value_type c=0.0,s=0.0;
		decode_givens(rho, c,s);
		std::size_t max_iter = std::min(std::max(i,k)+1,M.size2());
		for(std::size_t j = 0; j < max_iter; ++j){
			const double t1 = M(j,i);
			const double t2 = M(j,k);
			M(j,i) = c*t1 - s*t2;
			M(j,k) = s*t1 + c*t2;
		}
	}

	template< typename Matrix>
	void qr_iteration( Matrix & H, double tol=1e-16){
		#ifdef QR_ITERATION_OUTPUT
		std::cout << "H = " << print_matrix( H) << std::endl;
		#endif	
		std::size_t n = H.size1();
		typedef typename Matrix::value_type Value;
		typedef typename ublas::diagonal_adaptor< Matrix> 
							Diagonal_adapter;
		typedef typename ublas::scalar_matrix< Value> Scalar;
		typedef ublas::diagonal_adaptor< Scalar> Diagonal_scalar;
		std::vector< typename Matrix::value_type> givens(n-1, 0.0); 
		do{
			const Value a = H(n-2,n-2)+H(n-1,n-1);
			const Value b = std::sqrt(4*H(n-1,n-2)*H(n-2,n-1) + 
					std::pow((H(n-2,n-2)-H(n-1,n-1)),2));
			const Value e1 = (a+b)/2, e2 = (a-b)/2;
			const Value shift = (std::fabs(e1- H(n-1,n-1)) < 
					std::fabs(e2 - H(n-1,n-1)))? e1 : e2;
			
			Scalar mu_s(n,n,shift);
			const Diagonal_scalar mu( mu_s );
			Diagonal_adapter D(H);
			D -= mu;
			for (std::size_t i = 0; i < n-1; ++i){ 
				givens[i] = apply_givens_left(H,i,i+1); 
			}
			
			#ifdef DEBUG_QR_ITERATION
				std::cout << "Left Apply: H = " 
					<< print_matrix( H) << std::endl;
			#endif //DEBUG_QR_ITERATION
			
			for (std::size_t i = 0; i < n-1; ++i){ 
				apply_givens_right(H,givens[i],i,i+1); 	
			}
			
			#ifdef DEBUG_QR_ITERATION
				std::cout << "Right Apply H = " 
					<< print_matrix( H) << std::endl;
			#endif //DEBUG_QR_ITERATION
			
			//Unshift
			D += mu;
			
			#ifdef QR_ITERATION_OUTPUT
			std::cout.precision( 7);
			std::cout.setf( std::ios::fixed, std:: ios::floatfield);
			std::cout << "---- Shift: " << std::setw( 10) << shift  
			<< std::endl;
			std::cout << "---- Error: " << std::setw( 10)  
			<< std::fabs(H(n-1,n-2)) << std::endl;
			#endif //DEBUG_QR_ITERATION
		} while( std::fabs(H(n-1,n-2)) > tol);
		
	}
	template< typename Communicator, typename String>
	void debug_comm( const Communicator & c, const String & s){
		if( !c) { std::cerr << s << std::endl; }
	}
	template< typename Vector, typename Communicator>
	typename Vector::value_type compute_householder_vector( Vector & v, 
					 const Communicator & column_comm){
		typedef typename Vector::value_type Value;
		Value x = v[0];
		mpi::broadcast( column_comm, x, 0);
		Value beta = 1.0;
		const Value inner_prod = ublas::inner_prod(v,v);
		Value sigma = 0;
		mpi::all_reduce(column_comm, inner_prod, 
				sigma, std::plus< Value>()); 
		sigma -= x*x;
		if (sigma != 0){
			const Value mu = std::sqrt(x*x + sigma);
			(x <= 0)? x -= mu: x = -sigma/(x + mu);
			Value y = x*x;
			beta = (2.0*y)/(y+sigma);
			v /= x;
		}
		if( column_comm.rank()==0){ v[0]=1.0; }
		return beta;
	}
	
	template< typename Vector, typename Matrix, typename Communicator>
	void apply_householder_left( const typename Vector::value_type & beta, 
				     const Vector & v, 
				     Matrix & M, const std::size_t k, 
				     const Communicator & comm){
		const bool root = (comm.rank()==0);
	        if (beta != 0){
			std::cout << "r1: " << (k+1)*root << " -- " << M.size1() << std::endl
				  << "r2: " << k << " -- " << M.size2() << std::endl << std::flush; 
			ublas::range r1((k+1)*root, M.size1());
			ublas::range r2(k,  M.size2());
			ublas::matrix_range< Matrix> S(M, r1, r2);
			Vector w = ublas::prod< Vector>( ublas::trans( v), S);
			mpi::all_reduce( comm, w, w, std::plus< Vector>()); 
			S -= beta*ublas::outer_prod( v, w);
		}
	}

	template< typename Vector, typename Matrix, typename Communicator>
	void apply_householder_right( const typename Vector::value_type & beta, 
				      const Vector & v, 
				      Matrix & M, const std::size_t k, 
				      const Communicator & comm){ 
		const bool root = (comm.rank()==0);
		//if (beta != 0){
		//	ublas::range r1(0, M.size1());
		//	ublas::range r2(k*root, M.size2());
		//	ublas::matrix_range< Matrix> S(M, r1,r2);
		//	Vector w = ublas::prod<Vector>(S,v);
		//	mpi::all_reduce( comm, w, w, std::plus< Vector>()); 
		//	S -= beta*ublas::outer_prod( w, v);
		//}
	}
	template< typename Matrix_data>
	void apply_householder( Matrix_data & data, const std::size_t & k){ 
	typedef typename Matrix_data::Communicator Communicator; 
	typedef typename Matrix_data::Matrix Matrix; 
	typedef typename Matrix::value_type Value;
	typedef typename ublas::vector< Value> Vector;
	typedef typename std::pair< Value, Vector> Pair;
	const std::size_t block_col = block_column_index( k, data.row_length, 
								      data.n);
	const std::size_t column_index = local_column_index( k, data.row_length,
								 data.n);
	//const std::size_t row_index = k - data.first_row;
	const std::size_t col_root = data.block_col - block_col;
	const bool last_row = (k == data.last_row-1);
	Matrix & M = data.M;
	Vector vb_left( M.size1()+1, 0);
	if( !last_row ) { 
		mpi::broadcast( data.row_comm[ block_col], vb_left, 0); 
	}
        
	Vector vb_right( vb_left);
	vb_right.resize( M.size1()+1, M.size1()==M.size2());
	mpi::broadcast( data.col_comm[ block_col],  vb_right, col_root);

	Vector v_left = ublas::subrange(vb_left, 1, vb_left.size());
	Vector v_right = ublas::subrange(vb_right, 1, vb_right.size());
	const Value beta = vb_right[ 0];
	apply_householder_right( beta, v_right, M, 
				 column_index, data.col_comm[ block_col]);
	if ( !last_row ) { 
            apply_householder_left( beta, v_left, M, 
				    column_index, data.row_comm[ block_col]);
	}
	}

	template< typename Matrix_data>
	void dist_reduce_column( Matrix_data & data, const std::size_t & k){
	     typedef typename Matrix_data::Matrix Matrix;
	     typedef typename Matrix::value_type Value;
	     typedef typename ublas::matrix_column< Matrix> Matrix_column;
	     typedef typename ublas::vector< Value> Vector;
	     typedef ublas::slice Slice;
	     Matrix & M = data.M;
	     const bool on_last_col = (k == data.last_col-1);
	     const bool penult_blck_col = (data.block_col == data.row_length-2);
	     const bool serial = (on_last_col && penult_blck_col);
	     const std::size_t column_index = k - data.first_col; 
	     const std::size_t column_offset = data.diag()*(column_index+1);
	     Matrix_column col(M, column_index);
	     Vector v( col.size()-column_offset);
	     std::copy( col.begin()+column_offset, col.end(), v.begin());
	     const std::size_t comm_idx = 
	     block_column_index( k, data.row_length, data.n);
	     Value beta;
	     if ( serial ){
		serial::compute_householder_vector( v);
		beta = v[0]; v[0]=1;
	     } else {
	      const std::size_t col_comm_idx = comm_idx+on_last_col;
	      beta = compute_householder_vector( v, 
					data.col_comm[ col_comm_idx]);
	     }
	     std::cout << "k: " << k 
		       <<  " id: " << data.world.rank() 
		       << " v: " << v << std::endl;
	 
	     Vector w( v.size()+1);
	     std::copy( v.begin(), v.end(), w.begin()+1);
	     w[0] = beta;
	     mpi::broadcast( data.row_comm[ comm_idx], w, 0);
	     if (!serial){ 
	     mpi::broadcast( data.col_comm[ comm_idx], w, 0);
	     Vector hv_right = ublas::subrange(w, 1, w.size());
             apply_householder_right( beta, hv_right, M, column_index, 
				        data.col_comm[ comm_idx]);
	     }
	     apply_householder_left( beta, v, M, column_index, 
					data.row_comm[ comm_idx]);
	}
	template< typename Matrix_data>
	void hessenberg( Matrix_data & data){
	    const std::size_t id = data.world.rank();
	    const std::size_t n = data.n;
	    const std::size_t p = data.row_length;
	    //Algorithm 7.4.2 GVL
	    for (std::size_t k = 0; k < n-2; ++k){
		const std::size_t col_idx = block_column_index(k, p, n);
		if (data.block_row < col_idx || 
			data.block_col < col_idx){ return; }
		if (data.block_col == col_idx){
		     if(col_idx == data.row_length-1){
		     	serial::hessenberg( data.M);
			return;
		     }
		     if(data.diag() && k == data.last_col-1){ return; }
		     dist_reduce_column( data, k);  
		}else{ apply_householder( data, k); }
	   }
	}

	template< typename Matrix_data>
	void qr( Matrix_data & data){
			hessenberg( data);
			/*for (std::size_t i = M.size1(); i > 1; --i){
				Matrix_range R(M, Range (0, i), Range (0, i));
				qr_iteration( R);
			}*/
	}

} //end namespace parallel
} //end namespace t10
#endif //QR_DISTRIBUTED_H
