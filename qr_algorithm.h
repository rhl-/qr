#ifndef QR_ALGORITHM_H
#define QR_ALGORITHM_H
//#define DEBUG_QR_ITERATION //comment me out to turn off debug mode
#include <boost/timer.hpp>
#include <boost/numeric/ublas/banded.hpp> //for banded_adaptor
#include <boost/numeric/ublas/matrix.hpp> //for slices
#include <boost/numeric/ublas/matrix_proxy.hpp> //for slices
#include <boost/numeric/ublas/vector_proxy.hpp> //for slices
#include <boost/numeric/ublas/io.hpp> //io
#ifdef DEBUG_QR_ITERATION
#include <sstream>
#endif //endif DEBUG_QR_ITERATION
namespace ublas = boost::numeric::ublas;

namespace t10 {
	template< typename Matrix>
	std::string print_matrix( const Matrix & M){
		std::string str;
		std::stringstream ss;
		ss << M;
		str = ss.str();
		std::size_t lefts = str.find(std::string("["));
		std::size_t rights = str.find(std::string("]"));
		str.replace(lefts,rights-lefts+1,"");
		std::size_t found = str.find(std::string("(("));
		str.replace(found,2,std::string("["));
		found = str.find(std::string("))"));
		str.replace(found,2,std::string("]"));
		found = str.find(std::string("),"));
		while( found != std::string::npos){
			str.replace(found,2,std::string(";"));
			found = str.find(std::string("),"));
		}

		found = str.find(std::string("("));
		while( found != std::string::npos){
			str.replace(found,1,std::string(""));
			found = str.find(std::string("("));
		}


		return str;	
	}
	//GVL Section 5.1.9
	template< typename Value>
	void compute_givens(const Value & a, const Value & b, Value & c, Value & s){
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
		if (std::fabs(rho) < 1) {s = 2*rho; c = std::sqrt(1 -std::pow(s,2)); return;}
		c = 2/rho; s = std::sqrt(1-std::pow(c,2));
	} 
	
	template< typename Matrix>
	typename Matrix::value_type apply_givens_left( Matrix & M, 
						       const std::size_t i, const std::size_t k){
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
	void apply_givens_right( Matrix & M, const typename Matrix::value_type rho, 
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
		std::size_t n = H.size1();
		typedef typename ublas::diagonal_adaptor< Matrix> Diagonal_adapter;
		typedef typename ublas::scalar_matrix< typename Matrix::value_type> Scalar;
		typedef ublas::diagonal_adaptor< Scalar> Diagonal_scalar;
		std::vector< typename Matrix::value_type> givens(n-1, 0.0); 
		int counter=0;
		do{
		
		#ifdef DEBUG_QR_ITERATION
		std::cout << "---- Shift: " << H(n-1,n-1) << std::endl;
		#endif //DEBUG_QR_ITERATION
		Scalar mu_s(n,n,H(n-1,n-1));
		const Diagonal_scalar mu( mu_s );
		Diagonal_adapter D(H);
		D -= mu;
		for (std::size_t i = 0; i < n-1; ++i){ givens[i]=apply_givens_left(H,i,i+1); }
		#ifdef DEBUG_QR_ITERATION
		std::cout << "Left Apply: H = " << print_matrix( H) << std::endl;
		#endif //DEBUG_QR_ITERATION
		for (std::size_t i = 0; i < n-1; ++i){ apply_givens_right(H,givens[i],i,i+1); }
		#ifdef DEBUG_QR_ITERATION
		std::cout << "Right Apply H = " << print_matrix( H) << std::endl;
		#endif //DEBUG_QR_ITERATION
		//Unshift
		D += mu;
		#ifdef DEBUG_QR_ITERATION
		std::cout << "Shift: H = " << print_matrix( H) << std::endl;
		#endif //DEBUG_QR_ITERATION
		counter++;
		} while( counter < 100);
		
	}

	template< typename Vector>
	void compute_householder_vector( Vector & v){
		typedef typename Vector::value_type Value;
		Value x = v[0];
		Value beta = 1.0;
		const Value sigma = ublas::inner_prod(v,v) - x*x;
		if (sigma != 0){
			const Value mu = std::sqrt(x*x + sigma);
			(x <= 0)? x -= mu: x = -sigma/(x + mu);
			Value y = x*x;
			beta = (2.0*y)/(y+sigma);
			v /= x;
		}
		v[0] = beta;
	}
	
	template< typename Vector, typename Matrix>
	void apply_householder_left( const typename Vector::value_type & beta, 
				     const Vector & v, Matrix & M, const std::size_t k = 0){
		if (beta != 0){
			ublas::range r1(k+1, M.size1());
			ublas::range r2(k, M.size2());
			ublas::matrix_range< Matrix> S(M, r1,r2);
			S -= beta*ublas::outer_prod( v, ublas::prod<Vector>(ublas::trans(v),S));
		}
	}

	template< typename Vector, typename Matrix>
	void apply_householder_right( const typename Vector::value_type & beta, 
				      const Vector & v, Matrix & M, const std::size_t k = 0){
		if (beta != 0){
			ublas::range r1(0, M.size1());
			ublas::range r2(k+1, M.size2());
			ublas::matrix_range< Matrix> S(M, r1,r2);
			S -= beta*ublas::outer_prod( ublas::prod<Vector>(S,v), v);
		}
	}

	
	template< typename Matrix>
	void hessenberg( Matrix & M){
		typedef typename Matrix::value_type Value;
		typedef typename ublas::matrix_column< Matrix> Matrix_column;
		typedef typename ublas::vector< Value> Vector;
		typedef typename ublas::vector_range< Vector> Vector_range;
		const std::size_t n = M.size1();
		Value beta=0;
		//Algorithm 7.4.2 GVL
		for (std::size_t k = 0; k < n-2; ++k){
			//create a copy of the correct piece of the k^th column
			Vector vs = ublas::subrange(Matrix_column(M,k), k+1, n);
			//compute householder vector
			compute_householder_vector(vs);
			//hackery to store beta without extra space
			beta=vs[0]; vs[0]=1;
			apply_householder_left( beta, vs, M, k);
			apply_householder_right( beta, vs, M, k);
			//the algorithm hessenberg storing Q's would store v here
			//(betas would become a vector)
		}
	}

	template< typename Matrix>
	void qr( Matrix & M){
		typedef typename Matrix::value_type Value;
		switch(M.size1()) {
		case 0:
		case 1:
			return;
		case 2:
			const Value a = M(0,0)+M(1,1);
			const Value b = std::sqrt(4*M(1,0)*M(0,1) + std::pow((M(0,0)-M(1,1)),2));
			M(0,0) = (a + b)/2;
			M(1,1) = (a - b)/2;
			M(1,0) = 0;
			M(0,1) = 0;
			return;
		default:
			hessenberg(M);
			qr_iteration( M);
		}
	}

} //end namespace t10
#endif //QR_ALGORITHM_H
