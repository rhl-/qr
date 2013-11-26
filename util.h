#ifndef T10_UTIL_H                                          
#define T10_UTIL_H                                          

//BOOST MPI
#include <boost/mpi/group.hpp>

//Boost
#include <boost/unordered_map.hpp>

//STL
#include <algorithm>
#include <numeric>

namespace mpi = boost::mpi;

template< typename Stream, typename T>
Stream& operator<<( Stream & out, const  std::vector< T> & v){
	typedef typename std::vector< T> Vector;
	typedef typename Vector::const_iterator Iterator;
	for(Iterator i = v.begin(); i != v.end(); ++i){
		out << *i;
		if (i+1 != v.end()){ out << ", ";}
	}
	return out;
}


namespace t10 {

	template <class ForwardIterator, class T>
	void iota (ForwardIterator first, ForwardIterator last, T val){
	  while (first!=last) {
	    *first = val;
	    ++first;
	    ++val;
	  }
	}

	template< typename T>	
        std::pair<T, T> id_to_index(const T & proc_id, const T & p) {
                return std::make_pair(proc_id / p, proc_id % p);
        }

	template< typename T>	
        T index_to_id(const T & block_row, const T & block_col, const T & p) {
                return block_row * p + block_col;
        }

	template< typename Vector>
	void create_panel_vector( Vector & panel, const std::size_t k, 
						  const std::size_t row_length){
		typedef typename Vector::iterator Iterator;
		Iterator end = panel.begin()+(row_length-k);
		std::size_t a = k, b = k;
		for (Iterator i = panel.begin(); 
			      i != end; ++i,++a){
			*i = t10::index_to_id( a, k, row_length);
		}
	
		for(Iterator i = end; i != panel.end(); ++i){
		        *i = t10::index_to_id( k, ++b, row_length);
		}
	}

	template< typename Communicator, typename Vector> 
	boost::mpi::group create_group( const Communicator & world, 
						   const Vector & v){ 
		return world.group().include(v.begin(), v.end()); 
	}
	template< typename T>
	T row_id( const T & proc_id, const T & p){ return proc_id / p; }
	
	template< typename T>
	T col_id( const T & proc_id, const T & p){ return proc_id % p; }

	template< typename Matrix_data>
	void construct_communicators( Matrix_data & data){
		typedef typename Matrix_data::Communicator Communicator;
		typedef typename std::vector< std::size_t> Vector;
		typedef typename Vector::iterator Iterator;

		typedef typename std::vector< Communicator> Vector_comm;	
		typedef typename Vector_comm::iterator Comm_iterator;

		const Communicator & world = data.world;

		const std::size_t id = world.rank();
		const std::size_t world_size = world.size();
		const std::size_t row_length = std::sqrt(world_size);
		const std::size_t r = row_id( id, row_length);
		const std::size_t c = col_id( id, row_length);
		
		data.partner = index_to_id( c,r,  row_length);
		data.row_length = row_length;
		
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
		
		//communicators along row of this processor 
		//(for propagating right multiplication)
		Vector_comm  & row_comm = data.row_comm;
		//communicators along col of this processor 
		//(for propagating left multiplication)
		Vector_comm  & col_comm = data.col_comm;

		//now moving down the "diagonal" remove stuff from each
		//group making communicators	
		const std::size_t one = 1;
		for (std::size_t k = 0; k < std::max(data.block_col,one); ++k) {
		    Vector srow( row.begin()+k, row.end());
		    mpi::group row_group = t10::create_group(world, srow);
		    row_comm.push_back( mpi::communicator( world, row_group));
		    data.row.push_back( srow);
		}
 		
		for (std::size_t k = 0; k < std::max(data.block_row,one); ++k) {
		    Vector scol( col.begin()+k, col.end()); 
		    mpi::group col_group = t10::create_group(world, scol);
		    col_comm.push_back( mpi::communicator( world, col_group));
		    data.col.push_back( scol);
		}
	}

	std::size_t block_column_index(std::size_t k, std::size_t p, 
				  std::size_t number_of_rows) {
		
		const std::size_t avg_block_size = number_of_rows / p;
		// number of entries left over as before, which is 
		//number of blocks of lager size
		const std::size_t remaining = number_of_rows % p;
		// total number of columns which are within the larger blocks
		const std::size_t total = (avg_block_size + 1)*remaining;
		// check if k is within the first few larger blocks
		//otherwise first offset it by that many (matrix) columns and 
		//see how many blocks past the first "remaining" blocks is this 
		//block
		if (k < total) { return k / (avg_block_size + 1);}
		return (k-total) / avg_block_size + remaining;
	}

	template<typename _Matrix, typename _Communicator>
	struct Matrix_data {
		typedef _Matrix Matrix;
		typedef _Communicator Communicator;
		typedef typename boost::unordered_map< std::size_t,
						       std::size_t> Map;
		typedef typename std::vector< std::size_t> Vector;
		typedef typename std::vector< _Communicator> Vector_comm;

		bool special() const { 
			return (block_col == row_length-2) && 
				block_row==row_length-1;
		} 
		bool above() const { return block_row < block_col;}
		bool below() const { return block_col < block_row;}
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
		_Communicator ant_comm;
		
		Vector_comm col_comm;
		Vector_comm row_comm;
		std::vector< Vector> col; //for debugging
		std::vector< Vector> row; //for debugging
	}; // struct Matrix_data
}//end namespace t10

#endif//T10_UTIL_H
