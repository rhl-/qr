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
	void 	apply_householder_left(
				     const typename Vector::value_type & beta, 
				     Matrix & M, const std::size_t k, 
				     const Vector & v_left, 
				     const Communicator & comm,
				     const std::size_t id,
				     const bool row_flag=false,
				     const bool col_flag=true){
	        if (beta != 0){
			ublas::range rows((k+1)*row_flag, M.size1());
			ublas::range cols(k*col_flag,  M.size2());
			ublas::matrix_range< Matrix> S(M, rows, cols);
			Vector w = ublas::prod< Vector>( ublas::trans( v_left), S);
			mpi::all_reduce( comm, w, w, std::plus< Vector>()); 
			S -= beta*ublas::outer_prod( v_left, ublas::trans(w));
			}
	}

	template< typename Vector, typename Matrix, typename Communicator>
	void apply_householder_right( const typename Vector::value_type & beta, 
				      Matrix & M, const std::size_t k, 
				      const Vector & v_right, 
				      const Communicator & comm,
				      const std::size_t id,
				      const bool flag=false){ 
		if (beta != 0){
			ublas::range rows(0, M.size1());
			ublas::range cols((k+1)*flag, M.size2());
			ublas::matrix_range< Matrix> S(M, rows, cols);
			Vector w = ublas::prod<Vector>(S, v_right);
			mpi::all_reduce( comm, w, w, std::plus< Vector>());
			S -= beta*ublas::outer_prod( w, v_right);
		}
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
	const std::size_t col_root = data.block_col;
	const bool last_row = (k >= data.last_row-1);
	Matrix & M = data.M;
	Vector vb_left( M.size1()+1, 0);
	if( !last_row ){mpi::broadcast( data.row_comm[ block_col], vb_left, 0);}
        
	Vector vb_right( vb_left);
	vb_right.resize( M.size1()+1, M.size1()==M.size2());
	mpi::broadcast( data.col_comm[ 0],  vb_right, col_root);

	Vector v_left = ublas::subrange(vb_left, 1, vb_left.size());
	Vector v_right = ublas::subrange(vb_right, 1, vb_right.size());
	const Value beta = vb_right[ 0];
	const bool on_the_left = (block_col == data.block_col);
	const bool on_boundary_column = (k == (data.first_col-1));
	
	const bool on_penultimate_column = ((data.row_length-block_col)==2);
	const bool on_last_column = ((data.row_length-block_col)==1);
	const bool serial = (on_boundary_column && on_penultimate_column) 
				|| on_last_column;
	if( !serial){
	    apply_householder_right( beta, 
	   			 M, column_index, 
	   			 v_right, 
	   			 data.row_comm[ block_col+on_boundary_column],
	   			 data.world.rank(),
				 on_the_left);
	} else { 
		serial::apply_householder_right( beta, v_right, M, 
						 column_index*(on_the_left), 
						 on_the_left);
	}
	if (!last_row && serial){ 
		serial::apply_householder_left( beta, v_left, M, 0, 0); 
		return;
	}
	if ( !last_row ) {
	    const bool on_the_top = (block_col == data.block_row); 
	    apply_householder_left( beta, M, column_index,
				    v_left, data.col_comm[ block_col],
	   			    data.world.rank(),
				    on_the_top,
				    block_col >= data.block_col);
	}
	}

	template< typename Matrix_data>
	void dist_reduce_column( Matrix_data & data, const std::size_t & k){
	     typedef typename Matrix_data::Matrix Matrix;
	     typedef typename Matrix::value_type Value;
	     typedef typename ublas::matrix_column< Matrix> Matrix_column;
	     typedef typename ublas::vector< Value> Vector;
	     Matrix & M = data.M;

	     const std::size_t column_index = k - data.first_col; 
	     const std::size_t column_offset = data.diag()*(column_index+1);

	     Matrix_column col(M, column_index);
	     Vector v_left( col.size()-column_offset);
	     std::copy( col.begin()+column_offset, col.end(), v_left.begin());
	     const std::size_t comm_idx=block_column_index( k, data.row_length, data.n);
	     const std::size_t col_root = data.block_col;
	     Value beta;

	     const bool on_last_col = (k == data.last_col-1);
	     const bool penult_blck_col = (data.block_col == data.row_length-2);
	     const bool last_block_col = (data.block_col == data.row_length-1);
	     if ( (on_last_col && penult_blck_col) || last_block_col){
		serial::compute_householder_vector( v_left);
		beta = v_left[0]; v_left[0]=1;
	     } else {
	      const std::size_t col_comm_idx = comm_idx+on_last_col;
	      beta = compute_householder_vector( v_left, 
					data.col_comm[ col_comm_idx]);
	     }
	     Vector w( v_left.size()+1, beta);
	     std::copy( v_left.begin(), v_left.end(), w.begin()+1);
	     if(!last_block_col){mpi::broadcast( data.row_comm[ comm_idx], w, 0);}
	     if(!(on_last_col && penult_blck_col)) {
		 mpi::broadcast( data.col_comm[ 0], w, col_root); 
	     }
	     Vector v_right = ublas::subrange(w, 1, w.size());
	     const bool on_the_left = (comm_idx == data.block_col); 
	     if( !last_block_col){
	     apply_householder_right( beta, M, column_index, v_right,
				      	 data.row_comm[ comm_idx], 
					 data.world.rank(), 
					 on_the_left);
	     const bool on_the_top = (comm_idx == data.block_row); 
	     apply_householder_left( beta,  M, column_index,  v_left,
				     data.col_comm[ 0], 
				     data.world.rank(), on_the_top);
	    }else{
	      serial::apply_householder_right( beta, v_left, M, column_index);
	      serial::apply_householder_left( beta, v_left, M, column_index,); 
	    }
	}
	template< typename Matrix_data>
	void hessenberg( Matrix_data & data){
	    const std::size_t id = data.world.rank();
	    const std::size_t n = data.n;
	    const std::size_t p = data.row_length;
	    //Algorithm 7.4.2 GVL
	    for (std::size_t k = 0; k < n-2; ++k){
		const std::size_t col_idx = block_column_index(k, p, n);
		if ( data.block_col < col_idx ){ return; }
		if (data.block_col == col_idx && data.block_row >= col_idx){
		   if(data.diag() && k == data.last_col-1){ return; }
		   dist_reduce_column( data, k);  
		}else{ apply_householder( data, k); }
		std::cerr << print_matrix( data.M) << std::endl; 
	   }
		std::cerr << print_matrix( data.M) << std::endl; 
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
