#ifndef T10_UTIL_H                                          
#define T10_UTIL_H                                          

//STL
#include <algorithm>
#include <numeric>

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

	template< typename T>	
        std::pair<T, T> id_to_index(const T & proc_id, const T & p) {
                return std::make_pair(proc_id / p, proc_id % p);
        }

	template< typename T>	
        T index_to_id(const T & block_row, const T & block_col, const T & p) {
                return block_row * p + block_col;
        }

	template< typename T>
	T row_id( const T & proc_id, const T & p){ return proc_id / p; }
	
	template< typename T>
	T col_id( const T & proc_id, const T & p){ return proc_id % p; }

	std::size_t local_column_index(std::size_t k, std::size_t p,
					std::size_t number_of_rows) {
		//either k is within the larger blocks or not
		//if k is within larger blocks, return k mod larger size
		//else subtract total number of columns from the larger blocks
		//from k and then return k mod average size
		const std::size_t avg_block_size = number_of_rows / p;
		const std::size_t remaining = number_of_rows % p;
                // total number of columns which are within the larger blocks
                const std::size_t total = (avg_block_size + 1)*remaining;
		if (k < total) { return k % (avg_block_size + 1);}
                return (k-total) % avg_block_size ;
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
}//end namespace t10

#endif//T10_UTIL_H
