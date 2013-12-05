#ifndef T10_COMMUNICATORS_H
#define T10_COMMUNICATORS_H
//BOOST MPI
#include <boost/mpi/group.hpp>
#include <boost/mpi/timer.hpp>

//STL
#include <algorithm>
#include <numeric>

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


template< typename Communicator, typename Vector> 
boost::mpi::group create_group( const Communicator & world, 
					   const Vector & v){ 
	return world.group().include(v.begin(), v.end()); 
}

template< typename Matrix_data>
void construct_communicators( Matrix_data & data){
	std::cout << std::endl;
	mpi::timer t;
	typedef typename Matrix_data::Communicator Communicator;
	typedef typename Matrix_data::Vector_comm Vector_comm;
	typedef typename Vector_comm::iterator Comm_iterator;
	typedef typename std::vector< std::size_t> Vector;
	typedef typename Vector::iterator Iterator;

	const Communicator & world = data.world;
	const std::size_t id = world.rank();
	const std::size_t world_size = world.size();
	const std::size_t row_length = std::sqrt(world_size);
	const std::size_t r = row_id( id, row_length);
	const std::size_t c = col_id( id, row_length);
	
	data.partner = index_to_id( c,r,  row_length);
	data.row_length = row_length;
	std::cout  << "elapsed " << t.elapsed() << std::endl;
	t.restart();
	
	Vector row( row_length, id);
	Vector col( row_length, id);
	//create row and column mpi group indices
	//for proc and its partner
	for(Iterator i = row.begin(), j = col.begin(); 
		     i != row.end(); ++i, ++j){
		const std::size_t idx = std::distance(row.begin(), i);
		*i = r*row_length + idx; //row
		*j = idx*row_length + c; //col 
	}
	std::sort (row.begin(),   row.end()); 
	std::sort (col.begin(),   col.end()); 
	
	std::cout << "row: " << row << std::endl;
	std::cout << "col: " << col << std::endl;

	//communicators along row of this processor 
	//(for propagating right multiplication)
	Vector_comm  & row_comm = data.row_comm;
	//communicators along col of this processor 
	//(for propagating left multiplication)
	Vector_comm  & col_comm = data.col_comm;

	std::cout  << "elapsed " << t.elapsed() << std::endl;
	t.restart();
	//now moving down the "diagonal" remove stuff from each
	//group making communicators	
	const std::size_t one = 1;
	const std::size_t num_col_comm = data.block_row+1-
					(data.block_row == row_length-1);
	const std::size_t num_row_comm = data.block_col+1-
					(data.block_col == row_length-1);
	row_comm.reserve( num_row_comm);
	col_comm.reserve( num_col_comm); 
	for (std::size_t k = 0; k < num_row_comm; ++k) {
	    Vector srow( row.begin()+k, row.end());
	    std::cout << "srow: " << srow << std::flush;
	    mpi::group row_group = t10::create_group(world, srow);
	    std::cout << " group created" << std::endl;
	    std::cout << "col_comm size: " << row_comm.size() << "/" 
		      << num_row_comm << std::endl;
	    std::cout << " creating comm" << std::flush;
	    row_comm.emplace_back( world, row_group);
	    std::cout << " ...created!" << std::endl;
	}

	std::cout  << "elapsed " << t.elapsed() << std::endl;
	t.restart();
	for (std::size_t k = 0; k < num_col_comm; ++k) {
	    Vector scol( col.begin()+k, col.end()); 
	    std::cout << "scol: " << scol << std::flush;
	    mpi::group col_group = t10::create_group(world, scol);
	    std::cout << " group created" << std::endl;
	    std::cout << "col_comm size: " << col_comm.size() 
		      << "/" << num_col_comm << std::endl;
	    std::cout << " creating comm" << std::flush;
	    col_comm.emplace_back( world, col_group);
	    std::cout << " ...created!" << std::endl;
	}
	std::cout  << "elapsed " << t.elapsed() << std::endl;
	t.restart();
}
} //end namespace 10
#endif
