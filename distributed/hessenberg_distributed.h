#ifndef HESSENBERG_DISTRIBUTED_H
#define HESSENBERG_DISTRIBUTED_H

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
	    #ifdef HESSENBERG_OUTPUT
	std::cout << "apply_householder" << std::endl;
	#endif
	typedef typename Matrix_data::Communicator Communicator; 
	typedef typename Matrix_data::Matrix Matrix; 
	typedef typename Matrix::value_type Value;
	typedef typename ublas::vector< Value> Vector;
	typedef typename std::pair< Value, Vector> Pair;
	const std::size_t block_col = block_column_index( k, data.row_length, 
								      data.n);
	const std::size_t next_block_col=block_column_index( k+1, 
							     data.row_length, 
								      data.n);

	const std::size_t column_index = local_column_index( k, data.row_length,
								 data.n);
	const std::size_t col_root = data.block_col;
	Matrix & M = data.M;
	//Step 0: Check all the special cases
	const bool last_row = (k >= data.last_row-1);	
	const bool on_the_left = (block_col == data.block_col);
	const bool on_boundary_column = (block_col != next_block_col); 
	const bool on_penultimate_column = ((data.row_length-block_col)==2);
	const bool on_last_column = ((data.row_length-block_col)==1);
	const bool serial = (on_penultimate_column && on_boundary_column) || 
			    on_last_column;
	
	//Step 0.5: Create Left Householder Vector of appropriate size
	Vector vb_left( M.size1()+1, 0);

	//Step 1: If there is data to receive, grab it	
	if( !last_row ){ 
	    #ifdef HESSENBERG_OUTPUT
		std::cout << "row bcast" << std::endl;
		#endif
		mpi::broadcast( data.row_comm[ block_col], vb_left, 0);
	}
	Vector v_left = ublas::subrange(vb_left, 1, vb_left.size());
        
	//Step 1.5: create a Right Householder Vector of appropriate size
	Vector vb_right( vb_left);
	vb_right.resize( M.size1()+1, M.size1()==M.size2());

	//Step 2: Send/Receive the vector necessary for right applications
	    #ifdef HESSENBERG_OUTPUT
	std::cout << "col bcast" << std::flush;
	    #endif
	mpi::broadcast( data.col_comm[ 0],  vb_right, col_root);
	    #ifdef HESSENBERG_OUTPUT
	std::cout << "... done" << std::endl;
	    #endif
	Vector v_right = ublas::subrange(vb_right, 1, vb_right.size());
	const Value beta = vb_right[ 0];
	//Step 3: Apply the right householder matrix, sometimes in serial
	if( serial){
	    #ifdef HESSENBERG_OUTPUT
	     std::cout << "id: " << data.world.rank() 
	      	       << " applying: " << v_right << " on the right in serial"
		       << std::endl;
	    #endif
	     serial::apply_householder_right( beta, v_right, M, 
					     column_index*(on_the_left), 
				             on_the_left);
	} else {
	    const std::size_t comm_idx = block_col + on_boundary_column;
	    #ifdef HESSENBERG_OUTPUT
	    std::cout << "before right apply_householder" << std::endl; 
	    std::cout << "data.row_comm.size(): " << data.row_comm.size() << std::endl;
	    std::cout << "comm_idx: " << comm_idx << std::endl;
	   #endif 
	   //When we are on a boundary column one less processor is involved.
	    apply_householder_right( beta, 
	   			 M, column_index, 
	   			 v_right, 
	   			 data.row_comm[ comm_idx],
	   			 data.world.rank(),
				 on_the_left);
	}

	//Step 4: Apply the left householder matrix, sometimes in serial	
	if (!last_row && serial){ 
	    #ifdef HESSENBERG_OUTPUT
	      std::cout << "serial left apply_householder" << std::endl; 
             #endif
		serial::apply_householder_left( beta, v_left, M, 0, 0); 
		return;
	}
	if ( !last_row ) {
	    #ifdef HESSENBERG_OUTPUT
	    std::cout << "before left apply_householder" << std::endl; 
	    #endif
	    const bool on_the_top = (block_col == data.block_row); 
	    const std::size_t comm_idx = block_col + on_boundary_column;
	    apply_householder_left( beta, M, column_index,
				    v_left, data.col_comm[ comm_idx],
	   			    data.world.rank(),
				    on_the_top,
				    block_col >= data.block_col);
	}
	}
	
	template< typename Matrix_data>
	void dist_reduce_column( Matrix_data & data, const std::size_t & k){
	     #ifdef HESSENBERG_OUTPUT
	     std::cout << "dist_reduce_column" << std::endl;
	     #endif
	     typedef typename Matrix_data::Matrix Matrix;
	     typedef typename Matrix::value_type Value;
	     typedef typename ublas::matrix_column< Matrix> Matrix_column;
	     typedef typename ublas::vector< Value> Vector;
	     Matrix & M = data.M;

	     const std::size_t column_index = k - data.first_col; 
	     const std::size_t column_offset = data.diag()*(column_index+1);
	     const std::size_t comm_idx=
				block_column_index( k, data.row_length, data.n);
	     const std::size_t next_block_col=block_column_index( k+1, 
							     data.row_length, 
								      data.n);
	     const bool on_boundary_column = (comm_idx != next_block_col);
	     const std::size_t col_root = data.block_col;
	     //Step 0: Build Output Structures
	     Matrix_column col(M, column_index);
	     Vector v_left( col.size()-column_offset);
	     std::copy( col.begin()+column_offset, col.end(), v_left.begin());
	     Value beta;
	     //Step 0.5: Compute special cases
	     const bool on_last_col = (k == data.last_col-1);
	     const bool penult_blck_col = (data.block_col == data.row_length-2);
	     const bool last_block_col = (data.block_col == data.row_length-1);
	     const bool penult_boundary_col = (on_last_col && penult_blck_col);
	     const bool serial = penult_boundary_col || last_block_col;
	     //Step 1: Compute Householder Vector
	     //These special cases allow for serial computation
	     const std::size_t col_comm_idx = comm_idx+on_last_col;
	     if ( serial){
		serial::compute_householder_vector( v_left);
		beta = v_left[0]; v_left[0]=1;
	     } else { //Otherwise there is parallel computation
	      //However when on_last_col, one less processor is involved
	      beta = compute_householder_vector( v_left, 
					data.col_comm[ col_comm_idx]);
	     }
		#ifdef HESSENBERG_OUTPUT
	     std::cout << "id: " << data.world.rank() 
		       << " created: " 
		       << v_left << std::endl;	
	      #endif
	     //Step 2: Package and send householder vector 
	     //        used for left householder applications 
	     //        across current active row 
	     Vector w( v_left.size()+1, beta);
	     std::copy( v_left.begin(), v_left.end(), w.begin()+1);
	     if(!last_block_col){
		mpi::broadcast( data.row_comm[ comm_idx], w, 0);
	     }
	     //Step 3: Send/Receive Package containing
	     //        vector for right applications
	     if (!on_boundary_column){
		 mpi::broadcast( data.col_comm[ 0], w, col_root);
	     } 
	     Vector v_right = ublas::subrange(w, 1, w.size());
		#ifdef HESSENBERG_OUTPUT
	     std::cout << "id: " << data.world.rank() 
		       << " recv'd: " 
		       << v_right << std::endl;	
		#endif
	     //Step 4: Apply Householder on the right, sometimes in serial
	     const bool on_the_left = (comm_idx == data.block_col); 
	    if( !on_last_col && last_block_col){
		#ifdef HESSENBERG_OUTPUT
	     std::cout << "id: " << data.world.rank() 
	      	       << " applying: " << v_left << " on the right in serial"
		       << std::endl;
		#endif
		serial::apply_householder_right( beta, v_left, M, 
						 column_index, on_the_left);
	     }else if(!on_last_col && !penult_boundary_col){
		#ifdef HESSENBERG_OUTPUT
		std::cout << "id: " << data.world.rank() 
	      	      << " applying: " << v_right << " on the right in parallel"
		      << std::endl;
		#endif
	      apply_householder_right( beta, M, column_index, v_right,
				      	 data.row_comm[ comm_idx], 
					 data.world.rank(), 
					 on_the_left);
	     }
	     //Step 5: Apply householder on the left, sometimes in serial
	     if( serial) {
		if(penult_boundary_col){ 
			ublas::matrix_column< Matrix> pv(M, M.size2()-1);
			const Value sign = (pv[ 0] < 0)?-1:1;
			const Value beta = sign*ublas::norm_2(pv);
			pv( 0) = beta;
			std::fill(pv.begin()+1,pv.end(),0.0f);
			return; 
		}
		#ifdef HESSENBERG_OUTPUT
		std::cout << "id: " << data.world.rank() 
	      	       << " applying: " << v_left << " on the left in serial"
		       << std::endl;
		#endif
		serial::apply_householder_left( beta, v_left, M, column_index);
	     } else{
		#ifdef HESSENBERG_OUTPUT
		std::cout << "id: " << data.world.rank() 
	      	       << " applying: " << v_left << " on the left in parallel"
		       << std::endl;
		#endif
	     const bool on_the_top = (comm_idx == data.block_row);
	     apply_householder_left( beta,  M, column_index,  v_left,
				     data.col_comm[ col_comm_idx], 
				     data.world.rank(), on_the_top);
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
		#ifdef HESSENBERG_OUTPUT
		std::cout << "k: " << k << std::endl
			  << print_matrix( data.M) << std::endl; 
		#endif
		const bool on_last_col = (k == data.last_col-1); 
		if (data.block_col < col_idx){ return; }
		if (on_last_col && data.diag()) { return; }
		if (on_last_col && data.above() && col_idx < p-1) { return; } 
		if (data.block_col == col_idx && data.block_row >= col_idx){
		   dist_reduce_column( data, k);
	 	} else{ apply_householder( data, k); }
	   }
	}

} //end namespace parallel
} //end namespace t10
#endif //HESSENBERG_DISTRIBUTED_H
