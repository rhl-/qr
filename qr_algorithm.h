#ifndef QR_ALGORITHM_H
#define QR_ALGORITHM_H
#include <boost/timer.hpp>
#include <boost/numeric/ublas/matrix.hpp> //for slices
#include <boost/numeric/ublas/matrix_proxy.hpp> //for slices
#include <boost/numeric/ublas/vector_proxy.hpp> //for slices
#include <boost/numeric/ublas/io.hpp> //io
namespace ublas = boost::numeric::ublas;

namespace t10 {
	template< typename Matrix>
	void apply_givens_left( Matrix & M, const std::size_t i, const std::size_t j){}
	template< typename Matrix>
	void apply_givens_right( Matrix & M, const std::size_t i, const std::size_t j){}

	//TODO: Specialize for symmetric matrices
	template< typename Matrix, typename Matrix_list>
	void qr_iteration( Matrix & H, Matrix_list & Q){}
	
	template< typename Matrix>
	void qr_iteration( Matrix & H, double tol=1e-8){
		//TODO: Shift H	
		//TODO: Use tol to achieve deflation
		for (std::size_t i = 0; i < H.size1()-1; ++i){ apply_givens_left(H,i,i+1); }
		for (std::size_t i = 0; i < H.size1()-1; ++i){ apply_givens_right(H,i,i+1); }
		//TODO: Unshift H
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
