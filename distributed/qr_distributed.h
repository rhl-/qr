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

	template< typename Vector, typename Matrix_vector>
	void copy_vector( const Vector & source, Matrix_vector & dest){
	   std::copy( source.begin(), source.end(), dest.begin());
	}

	template< typename Matrix, typename Matrix_data>
	void recv_ghost_boundary_row( Matrix & H, 
					  Matrix_data & data, 
					  const std::size_t comm_idx){
	   typedef typename ublas::matrix_row< Matrix> Matrix_row;
	   const mpi::communicator & comm = data.acol_comm;
	   H.resize(H.size1()+1, H.size2(), true); 
	   Matrix_row last_row(H,H.size1()-1);
	   Vector _row;
	   comm.recv( comm.rank()+1, comm.rank(), _row);
	   copy_vector( _row, last_row);
	}

	template< typename Matrix, typename Matrix_data>
	void isend_ghost_boundary_row( Matrix & H, 
					  Matrix_data & data, 
					  const std::size_t comm_idx){
	   typedef typename ublas::matrix_row< Matrix> Matrix_row;
	   const mpi::communicator & comm = data.acol_comm;
	   Matrix_row first_row(H,0);
	   comm.isend( comm.rank()-1, comm.rank(), last_col);
	}

	template< typename Matrix, typename Matrix_data>
	void recv_ghost_boundary_column( Matrix & H, 
					  Matrix_data & data, 
					  const std::size_t comm_idx){
	   typedef typename ublas::matrix_column< Matrix> Matrix_column;
	   const mpi::communicator & comm = data.arow_comm[ comm_idx];
	   H.resize(H.size1(), H.size2()+1, true); 
	   Matrix_column last_col(H,H.size2()-1);
	   Vector _col;
	   comm.recv( comm.rank()+1, comm.rank(), _col);
	   copy_vector( _col, last_col);
	}

	template< typename Matrix, typename Matrix_data>
	void isend_ghost_boundary_column( Matrix & H, 
					  Matrix_data & data, 
					  const std::size_t comm_idx){
	        typedef typename ublas::matrix_column< Matrix> Matrix_column;
		const mpi::communicator & comm = data.arow_comm[ comm_idx];
		Matrix_column first_col(H,0);
		comm.isend( comm.rank()-1, comm.rank(), first_col);
	}

	template< typename Matrix, typename Matrix_data>
	void ghost_rows( Matrix & H, Matrix_data & data, 
			 const std::size_t comm_idx){
		if (data.block_row != 0){
		 isend_ghost_boundary_row( H, data, comm_idx);
		}

		if (!data.diag()){
		 recv_ghost_boundary_row( H, data, comm_idx);
		}
	}

	template< typename Matrix, typename Matrix_data>
	void ghost_columns( Matrix & H, Matrix_data & data, 
			    const std::size_t comm_idx,
			    const bool not_on_right){
		//Now send my first column to the left
		if(!data.diag()){
		 isend_ghost_boundary_column( H, data, comm_idx);
		}
		//if I have data to receive, grab it
		if(not_on_right){
		 recv_ghost_boundary_column( H, data, comm_idx);
		}
	}
	
	template< typename Matrix, typename Matrix_data, typename Givens>
	void diag_proc_qr( Matrix & H, Matrix_data & data, 
			   Givens & givens, const std::size_t comm_idx
			   const std::size_t block_index){
	 const mpi::communicator& comm = data.diag_comm[ comm_idx];
	 const bool not_on_right = (data.block_col != block_index);
	 
	 //wait until it is your turn to do left givens applies
	 if( comm.rank() ){ comm.recv( comm.rank()-1, mpi::any_tag); }
	 for (std::size_t i = 0; i < n-1; ++i){ 
	    givens[i]=apply_givens_left(H,i,i+1); 
	 }
	 //if you have someone to tell, bcast givens across row
	 if (not_on_right){ 
	  const mpi::communicator& row_comm = data.row_comm[comm_idx];
	  mpi::bcast( row_comm, givens, 0);
	 }
	 //TODO: isend message on diag_comm to start next row
	 //technically this is the RQ step, but, there is no harm in 
	 //doing this computation right now. so sans a column we can do this.
	 for (std::size_t i = 0; i < n-1; ++i){
	      //might need to ignore the last column
	      apply_givens_right(H,givens[i],i,i+1); 
	 }
	 //if you have someone to tell, bcast givens up column 
	 if(data.block_row != 0){
	  const mpi::communicator& col_comm = data.acol_comm;
	  mpi::bcast( col_comm, givens, col_comm.size()-1); 	
	 }
	}
	template< typename Matrix, typename Matrix_data, typename Givens>
	void diag_proc_rq( Matrix & H, Matrix_data & data, Givens & givens){
		//TODO: wait for column to arrive,
		//then do givens for last column
		//serial::apply_givens_right(H,givens[i],i,i+1); 
	}

	template< typename Matrix, typename Matrix_data, typename Givens>
	void offdiag_proc_rq( Matrix & H, Matrix_data & data, Givens & givens){
	  	const mpi::communicator& col_comm = data.acol_comm;
		//TODO: potentially check/wait for column to arrive
		//receive givens package column
		mpi::broadcast( col_comm, givens, col_comm.size()-1); 
		for (std::size_t i = 0; i < n-1; ++i){
		      serial::apply_givens_right(H,givens[i],i,i+1); 
		}
	}

	template< typename Matrix, typename Matrix_data, typename Givens>
	void offdiag_proc_qr( Matrix & H, Matrix_data & data, Givens & givens, 
			      const std::size_t comm_idx, 
			      const std::size_t block_index){
	  const mpi::communicator& row_comm = data.row_comm[comm_idx];
	  //Receive givens package from diagonal element to my left.
	  mpi::bcast( row_comm, givens, 0);
	  //TODO: write me --> apply_givens_left(H,givens);
	}
	template< typename Matrix, typename Matrix_data>
	void qr_iteration( Matrix & H, Matrix_data & data, 
			   const std::size_t block_index, 
			   const double tol=1e-16){
		typedef typename Matrix::value_type Value;
		#ifdef QR_ITERATION_OUTPUT
		std::cout << "H = " << t10::print_matrix( H) << std::endl;
		#endif
		const std::size_t n = H.size1();
		const bool not_on_right = (data.block_col != block_index);
		const bool not_on_bottom = data.block_row != block_index;
		
		std::vector< Value> givens(n-1, 0.0);
		const std::size_t comm_idx = data.row_length-1-block_index; 
		bool done = false;
		//do{

		  //Step 0: Compute Shift
		  if (data.diag()){ shift_matrix( data); }
		  
		  //Step 1: Ghost Rows
		  ghost_rows( H, data, comm_idx);

		  //Step 2: Do a dance
		  if(data.diag()) { 
		     diag_proc_qr( H, data, givens, comm_idx, block_index); 
		     #ifdef DEBUG_QR_ITERATION
		     std::cout << "Left Apply: H = " 
		     	  << t10::print_matrix( H) << std::endl;
		     #endif //DEBUG_QR_ITERATION
		     ghost_columns( H, data, comm_idx, not_on_right);
		     diag_proc_rq( H, data, givens); 
		  }
		  else {
		     offdiag_proc_qr( H, data, givens, comm_idx, block_index); 
		     #ifdef DEBUG_QR_ITERATION
		     std::cout << "Left Apply: H = " 
		     	  << t10::print_matrix( H) << std::endl;
		     #endif //DEBUG_QR_ITERATION
		     ghost_columns( H, data, comm_idx, not_on_right);
		     offdiag_proc_rq( H, data, givens); 
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
		  //} while( !done); 
		*/
	}
	
	template< typename Matrix_data>
	void qr( Matrix_data & data){
		typedef typename Matrix_data::Matrix Matrix;
		typedef typename ublas::matrix_range<Matrix> Matrix_range;
		typedef typename ublas::range Range;
		hessenberg( data);
		if( data.subdiag() || data.diag()){ 
			//TODO: send/recvone number
		}
		if( data.below()) { return; }
		const std::size_t begin = data.first_col;
		std::size_t end = data.last_col;
		for (std::size_t i = data.n; i > 1; --i){
			if (i < end){ end--; }
			if (end <= begin) { return; } 
			const std::size_t col_index = 
				block_column_index( i, data.row_length, data.n);
			Range r( begin, end);
			Matrix_range R(data.M, r, r);
			qr_iteration( R, data, col_index);
		}
	}

} //end namespace parallel
} //end namespace t10
#endif //QR_DISTRIBUTED_H
