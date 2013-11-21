#ifndef T10_UTIL_H                                          
#define T10_UTIL_H                                          
namespace t10 {
		template<typename _Matrix, typename Communicator>
	        struct Matrix_data {
			//Internal Matrix Type
	        	typedef typename _Matrix Matrix;
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
			
			/*
 			* Below is a list of communicators
 			*/ 	
			Communicator row_comm;
			Communicator col_comm;
			Communicator antidiag_comm;
		}; // struct Matrix_data

}//end namespace t10
#endif//T10_UTIL_H
