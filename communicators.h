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
	std::cout << std::endl;
	mpi::timer t;
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
	const std::size_t num_col_comm = data.block_row+1-
					(data.block_row == row_length-1);
	const std::size_t num_row_comm = data.block_col+1-
					(data.block_col == row_length-1);
	row_comm.reserve( num_row_comm);
	for (std::size_t i = 0; i < row_length; ++i){
		Vector row(row_length,0);
		std::iota( row.begin(), row.end(), i*row_length);
		for( std::size_t j = 0; j < row_length-1; ++j){
			Vector tmp( row.begin()+j, row.end());
			mpi::group row_group = t10::create_group( world,
								  row.begin()+j,
								  row.end());
			mpi::communicator comm( world, row_group);
			if( std::binary_search( row.begin()+j, row.end(), id)){ 
			   std::cout << "cur row: " << tmp << std::flush;
			   row_comm.emplace_back( comm, mpi::comm_attach);
			   std::cout << "... created" << std::endl;
			}
		}
	}
	col_comm.reserve( num_col_comm);
	for (std::size_t i = 0; i < row_length; ++i){
		Vector col(row_length,0);
		std::iota( col.begin(), col.end(), 0);
		t10::vec_multiply( col.begin(), col.end(), 3);
		t10::vec_add( col.begin(), col.end(), i);
		for( std::size_t j = 0; j < row_length-1; ++j){
			Vector tmp( col.begin()+j, col.end());
			mpi::group col_group = t10::create_group( world,
								  col.begin()+j,
								  col.end());
			mpi::communicator comm( world, col_group);
			if( std::binary_search( col.begin()+j, col.end(), id)){ 
			    std::cout << "cur col: " << tmp << std::flush;
			    col_comm.emplace_back( comm, mpi::comm_attach);
			    std::cout << "... created" << std::endl;
			}
		
		}
	}
}
} //end namespace 10
#endif
