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

	//by default these should be empty.	
	Vector_comm col_comm;
	Vector_comm row_comm;
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
	
	//communicators along row of this processor 
	//(for propagating right multiplication)
	Vector_comm  & row_comm = data.row_comm;
	//communicators along col of this processor 
	//(for propagating left multiplication)
	Vector_comm  & col_comm = data.col_comm;

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
		if (!row_comm[i-1]) { std::cerr << "row comm is invalid!" << std::endl;}
		mpi::communicator next = row_comm[ i-1].split( row_comm[i-1].rank()>0);
		if (next.size() == 1){ break; }
		row_comm.emplace_back( next, mpi::comm_duplicate);
	}
	for( std::size_t i = 1; i < row_length; ++i){
		if (!col_comm[ i-1]) { std::cerr << "col comm is invalid!" << std::endl;}
		mpi::communicator next = col_comm[ i-1].split( col_comm[i-1].rank() > 0);
		if (next.size() == 1){ break; }
		col_comm.emplace_back( next, mpi::comm_duplicate); 
	}
}
} //end namespace 10
#endif
