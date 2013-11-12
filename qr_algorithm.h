#ifndef QR_ALGORITHM_H
#define QR_ALGORITHM_H
#include <boost/timer.hpp>
#include <boost/numeric/ublas/banded.hpp> //for banded_adaptor
#include <boost/numeric/ublas/matrix.hpp> //for slices
#include <boost/numeric/ublas/matrix_proxy.hpp> //for slices
#include <boost/numeric/ublas/vector_proxy.hpp> //for slices
#include <boost/numeric/ublas/io.hpp> //io
#include <sstream>
namespace ublas = boost::numeric::ublas;

namespace t10 {
	template< typename Matrix, typename String>
	String print_matrix( const Matrix & M, String & str){
		std::stringstream ss;
		ss << M;
		str = ss.str();
		std::size_t lefts = str.find(std::string("["));
		std::size_t rights = str.find(std::string("]"));
		str.replace(lefts,rights-lefts+1,"");
		std::size_t found = str.find(std::string("(("));
		str.replace(found,2,std::string("\n("));
		found = str.find(std::string("))"));
		str.replace(found,2,std::string(")\n"));
		found = str.find(std::string("),"));
		while( found != std::string::npos){
			str.replace(found,2,std::string(")\n"));
			found = str.find(std::string("),"));
		}
		return str;	
	}
	//GVL Section 5.1.9
	template< typename Value>
	void compute_givens(const Value & a, const Value & b, Value & g0, Value & g1){
			if (b == 0) { g0 = 1; g1 = 0; return; }
			if( std::fabs(a) < std::fabs(b)){
				const Value rat = a/b;
				g1 = 1/ std::sqrt(1 + std::pow(rat, 2));
				g0 = -rat * g1;
				return;
			}
			const Value rat = b/a;
			g0 = 1/ std::sqrt(1 + std::pow(rat, 2));
			g0 = -rat * g0;
	}

	template< typename Matrix>
	void apply_givens_left( Matrix & M, const std::size_t i){
		typename Matrix::value_type g0=0.0,g1=0.0;
		compute_givens(M(i,i),M(i+1,i), g0,g1);
		for(std::size_t k = i; k < M.size2(); ++k){
			const double ti = M(i,k);
			const double tj = M(i+1,k);
			M(i,k) = g0*ti - g1*tj;
			M(i+1,k) = g1*ti + g0*tj;
		}
			
	}
	template< typename Matrix>
	void apply_givens_right( Matrix & M, const std::size_t i){
		typename Matrix::value_type g0=0.0,g1=0.0;
		compute_givens(M(i,i),M(i+1,i), g0,g1);
		for(std::size_t k = i; k < M.size1(); ++k){
			const double ti = M(k,i);
			const double tj = M(k,i+1);
			M(k,i) = g0*ti + g1*tj;
			M(k,i+1) = -g1*ti + g0*tj;
		}
	}

	//TODO: Specialize for symmetric matrices
	template< typename Matrix, typename Matrix_list>
	void qr_iteration( Matrix & H, Matrix_list & Q){}
	
	template< typename Matrix>
	void qr_iteration( Matrix & H, double tol=1e-16){
		std::size_t n = H.size1();
		typedef typename ublas::diagonal_adaptor< Matrix> Diagonal_adapter;
		typedef typename ublas::scalar_matrix< typename Matrix::value_type> Scalar;
		typedef ublas::diagonal_adaptor< Scalar> Diagonal_scalar;
		//Choose and apply shift
		std::string sout;
		do{
		std::cout << "----" << std::endl;
		Scalar mu_s(n,n,H(n-1,n-1));
		const Diagonal_scalar mu( mu_s );
		Diagonal_adapter D(H);
		D -= mu;
		//TODO: Use tol to achieve deflation
		for (std::size_t i = 0; i < n-1; ++i){ apply_givens_left(H,i); 
		}
		std::cout << "Left Apply: H = " << print_matrix( H, sout) << std::endl;
		for (std::size_t i = 0; i < n-1; ++i){ apply_givens_right(H,i); }
		std::cout << "Right Apply H = " << print_matrix( H, sout) << std::endl;
		//Unshift
		D += mu;
		std::cout << "Shift: H = " << print_matrix( H, sout) << std::endl;
		} while( std::fabs(H(n-1,n-2)) > tol);
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

	template< typename Matrix, typename Matrix_list>
	void hessenberg( Matrix & M, Matrix_list & Q){}
	
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
		//handle edge cases gracefully.
		if (M.size1() <= 2){return;}
			//TODO: compute answer instantly;

		hessenberg(M);
		qr_iteration( M);
		//TODO: extract D
	}

	template< typename Matrix>
	void qr( Matrix & M, Matrix & V){
		//handle edge cases gracefully.
		if (M.size1() <= 2){
			//TODO: compute answer instantly;
			return;
		}
		std::vector< Matrix> Q; 
		hessenberg(M, Q);
		qr_iteration( M, Q);
		//TODO: compute V
	}
} //end namespace t10
#endif //QR_ALGORITHM_H
