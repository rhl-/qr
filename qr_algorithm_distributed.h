#ifndef QR_ALGORITHM_H
#define QR_ALGORITHM_H
#define QR_HOUSE_DEBUG
#define QR_ITERATION_OUTPUT
//#define DEBUG_QR_ITERATION //comment me out to turn off debug mode
#include <boost/timer.hpp>
#include <boost/numeric/ublas/banded.hpp> //for banded_adaptor
#include <boost/numeric/ublas/matrix.hpp> //for slices
#include <boost/numeric/ublas/matrix_proxy.hpp> //for slices
#include <boost/numeric/ublas/vector_proxy.hpp> //for slices
#include "io.h"
#include "util.h"

#ifdef DEBUG_QR_ITERATION
#endif //endif DEBUG_QR_ITERATION
namespace ublas = boost::numeric::ublas;

namespace t10 {
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
	int sign( const T & x){return (x < 0) ? -1 : 1;}
	
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

	//TODO: Specialize for symmetric matrices
	
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
	//TODO: MPI group
	template< typename Vector, typename Communicator>
	typename Vector::value_type compute_householder_vector( Vector & v, 
					 const Communicator & column_comm){
		typedef typename Vector::value_type Value;
		Value x = v[0];
		#ifdef QR_HOUSE_DEBUG
		std::cout << "Proc: " << column_comm.rank()
			  << " Initial: " << v << std::endl;
		#endif //QR_HOUSE_DEBUG
		//root process should be zero, according to at least OpenMPI
		//documentation
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
		#ifdef QR_HOUSE_DEBUG
		std::cout << "Proc: " << column_comm.rank()
			  << " Beta: " << beta << " " << v << std::endl;
		#endif //QR_HOUSE_DEBUG
		return beta;
	}
	
	template< typename Vector, typename Matrix>
	void apply_householder_left( const typename Vector::value_type & beta, 
				     const Vector & v, Matrix & M ){
		/*
		if (beta != 0){
			Vector w = ublas::prod< Vector>(M,v);
			mpi::all_reduce( col_comm, w, w, std::plus< Vector>()); 
			M -= beta*ublas::outer_prod( v,w);
		}*/
	}

	template< typename Vector, typename Matrix>
	void apply_householder_right( const typename Vector::value_type & beta, 
				      const Vector & v, 
				      Matrix & M, const std::size_t k = 0){
		/*
		if (beta != 0){
			Vector w = ublas::prod<Vector>(M,v);
			mpi::all_reduce( row_comm, w, w, std::plus< Vector>()); 
			M -= beta*ublas::outer_prod( w, v);
		}*/
	}
	
	template< typename Matrix_data, typename T>
	std::pair< T, T> compute_roots(  const Matrix_data & data, 
					const T & block_col){
		const std::size_t lroot = index_to_id( data.block_row,block_col,
						       data.row_length);
		const std::size_t rroot = index_to_id( data.block_col,block_col,
						       data.row_length);
		return std::make_pair( lroot, rroot);
	}
	template<typename Map>
	std::string print_map(const Map & map){
		typedef typename Map::const_iterator Iterator;
		std::stringstream ss;
		ss << "{" << std::endl;
		for(Iterator i = map.begin(); i != map.end(); ++i){
			ss << i->first << " --> " << i->second << std::endl; 
		}
		ss << "}" << std::endl;
		return ss.str();
	}
	template< typename Matrix_data, typename T>
	bool get_right_rank( const Matrix_data & data, 
			    const T & block_col, std::size_t & rank){
		typedef typename Matrix_data::Map_vector::value_type Map;
		typedef typename Map::const_iterator Iterator;
		const std::size_t id = data.world.rank();
		const Map & map = data.right_comm_map[ block_col];
		std::cout << id << " rcms " << data.right_comm_map.size()
			  << " map_size: " << map.size()
			  << std::endl
			  << print_map( map);
		Iterator pair = map.find( data.world.rank());
		if (pair == map.end()){
			return false; 
		}
		rank = pair->second;
		return true; 
	}

	template< typename Matrix_data, typename T>
	T get_left_rank( const Matrix_data & data, const T & block_col){
		typedef typename Matrix_data::Map_vector::value_type Map;
		typedef typename Map::const_iterator Iterator;
		const std::size_t id = data.world.rank();
		const Map & map = data.left_comm_map[ block_col];
		Iterator pair = map.find( data.world.rank());
		if (pair == map.end()){ 
			std::cerr << id << " bug exists get_left_rank error" 
			          << std::endl;
			return 666; 
		}
		return pair->second; 
	}

	template< typename Matrix_data>
	void hessenberg( Matrix_data & data){
		typedef typename Matrix_data::Communicator Communicator;
		typedef typename Matrix_data::Matrix Matrix;
		typedef typename Matrix::value_type Value;
		typedef typename ublas::matrix_column< Matrix> Matrix_column;
		typedef typename ublas::vector< Value> Vector;
		typedef typename ublas::vector_range< Vector> Vector_range;
		typedef typename std::pair< std::size_t, std::size_t> Pair;
		Matrix & M = data.M;
		const std::size_t n = data.n;
		const std::size_t id = data.world.rank();
		const std::size_t blocks = std::sqrt( data.world.size());
		bool ready_to_load_balance=false;
		//Algorithm 7.4.2 GVL
		for (std::size_t k = 0; k < n-2; ++k){
			const std::size_t block_col = 
			t10::block_column_index( k, blocks, n);
			Communicator left_comm = data.left_comm[ block_col];
			Communicator right_comm = data.right_comm[ block_col];
			Pair root = compute_roots( data, block_col);
			if( k < data.first_col){
				Value beta=0.0;
				Vector v_left, v_right;
				//TODO: if possibly these should be Ibcast,
				//then wait_some() and as data comes in we call
				//the functions below
				mpi::broadcast(left_comm,  v_left, root.first);
				mpi::broadcast(right_comm, v_right,root.second);
				std::cout << id << " v_left: " 
					  << v_left << std::endl; 
				std::cout << id << " v_right: " 
					  << v_right << std::endl; 
				/*
				apply_householder_left( beta, 
							vs_left, M, row_comm);
				apply_householder_right( beta, 
							 vs_right, M, col_comm);
				*/
			}
			else if(data.below() && k < data.last_col-1){ 
				const std::size_t col_idx = k-data.first_col;
				Matrix_column col(M,col_idx);
				std::size_t offset = data.diag(); 
				Vector vs = ublas::subrange( col,
							     col_idx+offset, 
							     M.size1());
				const Communicator & cc = data.l_col_comm;
				Value beta = compute_householder_vector( vs,cc);
				//TODO: attach beta to vs
				std::size_t right_rank; 
				if(!get_right_rank( data, block_col, 
							  right_rank)){
					std::cout << "k = " << k << std::endl
						  << "id = " << id << std::endl
						  << "block_col = " << block_col
						  << std::endl;
					std::exit( -1);
				}
				mpi::broadcast( right_comm, vs, right_rank);
				std::cout << id << " just send hv" << std::endl;
			}
			//the last column of every block has a special case
			else if( data.below() && !data.diag() 
					&& k == data.last_col){
				const std::size_t col_idx = k-data.first_col;
		
				Matrix_column col(M,M.size2()-1);
				std::size_t offset =(data.s_col_comm.rank()==0);
				Vector vs = ublas::subrange( col,col_idx+offset,
							     M.size1());

				const Communicator & cc = data.s_col_comm;
				Value beta = compute_householder_vector( vs,cc);
				//TODO: attach beta to vs
				std::size_t right_rank; 
				if(!get_right_rank( data, block_col, 
							  right_rank)){
					std::cout << "k = " << k << std::endl
						  << "id = " << id << std::endl
						  << "block_col = " << block_col
						  << std::endl;
					std::exit( -1);
				}
				mpi::broadcast( right_comm, vs, right_rank);
				std::cout << id << " just send hv" << std::endl;
			}
			else if( !data.diag() ){
				ready_to_load_balance = true;
				break;
			}
			//JUST FOR DEBUGGING
			break;
		}
		if (ready_to_load_balance){
			if(data.above()){
				ublas::range r1(0,M.size1()/2);
				ublas::range r2(0, M.size2());
				ublas::matrix_range< Matrix> S(M, r1,r2);
				//TODO: Fix serialization issues
				//data.world.send(data.partner,0,S.data());
			}
			else{
				M.resize(M.size1(),M.size2()/2);
				//TODO: Fix serialization issues
				//data.world.recv(data.partner,0,M); 
			}
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

} //end namespace t10
#endif //QR_ALGORITHM_H
