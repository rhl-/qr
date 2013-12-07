#ifndef T10_COMMUNICATORS_H
#define T10_COMMUNICATORS_H
//BOOST MPI
#include <boost/mpi/group.hpp>
#include <boost/mpi/timer.hpp>

//STL
#include <algorithm>
#include <numeric>
#include <functional>   // std::multiplies

#include "util.h"

namespace mpi = boost::mpi;
namespace t10 {
template<typename _Matrix, typename _Communicator>
struct Matrix_data {
	typedef _Matrix Matrix;
	typedef _Communicator Communicator;
	typedef typename std::vector< std::size_t> Vector;
	typedef typename std::vector< _Communicator> Vector_comm;

	bool above() const { return block_row < block_col;  }
	bool below() const { return block_col < block_row;  }
	bool diag()  const { return block_row == block_col; }

	//Internal Matrix Type
	Matrix M;
	//Overall matrix size
	std::size_t n;

	//indices into the theoretical larger matrix
	std::size_t first_row;
	std::size_t last_row;
	std::size_t first_col;
	std::size_t last_col;
	std::size_t block_row;
	std::size_t block_col;
	std::size_t partner;
	std::size_t row_length;
	
	/*
	* Below is a list of communicators
	*/ 	
	_Communicator world;
	_Communicator above_col;
	_Communicator above_comm;

	//by default these should be empty.	
	Vector_comm col_comm;
	Vector_comm row_comm;
	Vector_comm diag_comm;
	Vector_comm arow_comm;
}; // struct Matrix_data


template< typename Communicator, typename Iterator> 
boost::mpi::group create_group( const Communicator & world, 
				const Iterator & begin,
				const Iterator & end){ 
	return world.group().include(begin, end); 
}

template< typename Iterator, typename T>
Iterator vec_add( Iterator begin, Iterator end, const T & t){
	return std::transform( begin, end, begin, 
			       std::bind2nd( std::plus<T>(), t) );
}


template< typename Iterator, typename T>
Iterator vec_multiply( Iterator begin, Iterator end, const T & t){
	return std::transform( begin, end, begin, 
			       std::bind2nd( std::multiplies<T>(), t) );
}



template< typename Matrix_data>
void construct_communicators( Matrix_data & data){
	typedef typename Matrix_data::Communicator Communicator;
	typedef typename Matrix_data::Vector_comm Vector_comm;
	typedef typename std::vector< std::size_t> Vector;

	const Communicator & world = data.world;
	const std::size_t id = world.rank();
	const std::size_t world_size = world.size();
	const std::size_t row_length = std::sqrt(world_size);
	const std::size_t r = row_id( id, row_length);
	const std::size_t c = col_id( id, row_length);
	
	data.partner = index_to_id( c, r,  row_length);
	data.row_length = row_length;
	
	Communicator & above_comm = data.above_comm;
	Communicator & above_col = data.above_col;
	//communicators along row of this processor 
	//(for propagating right multiplication)
	Vector_comm  & row_comm = data.row_comm;
	//communicators along col of this processor 
	//(for propagating left multiplication)
	Vector_comm  & col_comm = data.col_comm;
	Vector_comm & diag_comm = data.diag_comm;
	Vector_comm & arow_comm = data.arow_comm; 

	//Communicator k contains {1,2,..,p-k} diagonal processors

	//communicators for QR iteration
	mpi::communicator _abve_comm = world.split(data.above() || data.diag());
	mpi::communicator _diag_comm  = world.split( data.diag());
	mpi::communicator first_row = above_comm.split( data.block_row);
	mpi::communicator _col = above_comm.split( data.block_col);
	if( _col.size() > 1){ above_col = _col; } 
	if (data.diag()) {  
	      diag_comm.reserve(row_length-data.block_col);
	      diag_comm.emplace_back( _diag_comm, mpi::comm_duplicate);
	      for( std::size_t i = 1; i < row_length-1; ++i){
		const mpi::communicator & prev = diag_comm[i-1];
	        const bool key = (prev.rank() == (prev.size()-1));
	      	mpi::communicator next = diag_comm[ i-1].split( key);
	        if (key){ break; }
	      	diag_comm.emplace_back( next, mpi::comm_duplicate);
	      }
	}
	//now build arow_comm
	if ( first_row.size() > 1){
	   arow_comm.reserve(row_length-data.block_col);
	   arow_comm.emplace_back( first_row, mpi::comm_duplicate);
	   for( std::size_t i = 1; i < row_length-1; ++i){
		const mpi::communicator & prev = arow_comm[i-1];
		const bool key = (prev.rank() == prev.size()-1);
		mpi::communicator next = arow_comm[i-1].split(key);
		if(next.size() == 1) { break; }
		arow_comm.emplace_back( next, mpi::comm_duplicate);	
	   }
	}
	//now moving down the "diagonal" remove stuff from each
	//group making communicators	
	const std::size_t one = 1;
	row_comm.reserve( row_length-1);
	col_comm.reserve( row_length-1);
	mpi::communicator row = world.split( data.block_row);
	mpi::communicator col = world.split( data.block_col);
	row_comm.emplace_back( row, mpi::comm_duplicate);
	col_comm.emplace_back( col, mpi::comm_duplicate);
	std::array< std::size_t, 1> ids = {0};
	for( std::size_t i = 1; i < row_length; ++i){
		const bool key = row_comm[i-1].rank()>0;
		mpi::communicator next = row_comm[ i-1].split( key);
		if (next.size() == 1){ break; }
		row_comm.emplace_back( next, mpi::comm_duplicate);
	}
	for( std::size_t i = 1; i < row_length; ++i){
		const bool key = col_comm[i-1].rank()>0;
		mpi::communicator next = col_comm[ i-1].split( key);
		if (next.size() == 1){ break; }
		col_comm.emplace_back( next, mpi::comm_duplicate); 
	}
}
} //end namespace t10
#endif //T10_COMMUNICATORS_H
