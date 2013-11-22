#ifndef T10_UTIL_H                                          
#define T10_UTIL_H                                          
#define UTIL_DEBUG 
//Boost MPI
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
namespace t10 {

	template< typename T>
	T row_id( const T & proc_id, const T & p){ return proc_id / p; }
	
	template< typename T>
	T col_id( const T & proc_id, const T & p){ return proc_id % p; }
	
	template< typename T>	
        std::pair<T, T> id_to_index(const T & proc_id, const T & p) {
                return std::make_pair(proc_id / p, proc_id % p);
        }

	template< typename T>	
        T index_to_id(const T & block_row, const T & block_col, const T & p) {
                return block_row * p + block_col;
        }

	template< typename Matrix_data>
	void construct_communicators( Matrix_data & data){
		typedef typename Matrix_data::Communicator Communicator;
		typedef typename std::vector< std::size_t> Vector;
		typedef typename Vector::iterator Iterator;
		const Communicator & world = data.world;
		const std::size_t id = world.rank();
		const std::size_t world_size = world.size();
		const std::size_t row_length = std::sqrt(world_size);
		const std::size_t r = row_id( id, row_length);
		const std::size_t c = col_id( id, row_length);
		const std::size_t diag1 = std::min(r,c);
		const std::size_t diag2 = std::max(r,c);
		data.partner = index_to_id( c,r,  row_length);
		const std::size_t anti_m = std::min(r+c, row_length-1);
		const std::size_t anti_size = 2*anti_m-(r+c) + 1;
		Vector ant( anti_size, id);
		Vector row( row_length, id);
		Vector col( row_length, id);
		Vector pan( 2*(row_length-diag1)-1, id);

		//create panel indices
		Iterator end = pan.begin()+(row_length-diag1);
		std::size_t a = diag1;
		std::size_t b = diag1;
		for (Iterator i = pan.begin(); i != end; ++i,++a){
				*i = index_to_id(a,diag1, row_length);
		}

		for(Iterator i = end; i != pan.end(); ++i){
				*i = index_to_id(diag1,++b, row_length);
		}

		//create row and column mpi group indices
		Iterator j = col.begin();
		for(Iterator i = row.begin(); i != row.end(); ++i, ++j){
			const std::size_t k = std::distance(row.begin(), i);
			*i = r*row_length + k;
			*j = k*row_length + c;
		}

		//create antidiagonal group indices
		for (Iterator k = ant.begin(); k != ant.end(); ++k){
		       const std::size_t i = std::distance(ant.begin(), k);
			*k = index_to_id((r+c)-anti_m+i,anti_m-i, row_length);
		}
		#ifdef UTIL_DEBUG
			std::cout << "Processor: " << world.rank() <<  std::endl
				  << "Partner: " << data.partner << std::endl
				  << "row: " << row << std::endl 
				  << "col: " << col << std::endl 
				  << "ant: " << ant << std::endl
				  << "pan: " << pan << std::endl; 
		#endif
		mpi::group row_group = world.group().include( row.begin(), 
								row.end());
		mpi::group col_group = world.group().include( col.begin(), 
								col.end());
		mpi::group ant_group = world.group().include( ant.begin(), 
								ant.end());
		mpi::group pan_group = world.group().include( pan.begin(), 
								pan.end());

		data.row_comm = mpi::communicator(world, row_group);
		data.col_comm = mpi::communicator(world, col_group);
		data.ant_comm = mpi::communicator(world, ant_group);
		data.pan_comm = mpi::communicator(world, pan_group);
		#ifdef UTIL_DEBUG
		std::cout << "Processor: " << world.rank() << " is " 
			  << std::endl
			  << data.row_comm.rank() << " of " 
			  << data.row_comm.size() << " in row_comm" 
			  << std::endl
			  
			  << data.col_comm.rank() << " of " 
			  << data.col_comm.size() << " in col_comm" 
			  << std::endl
 			  
			  << data.ant_comm.rank() << " of " 
			  << data.ant_comm.size() << " in ant_comm" 
			  << std::endl

			  << data.pan_comm.rank() << " of " 
			  << data.pan_comm.size() << " in pan_comm" 
			  << std::endl
			  << "-------------------------" 
			  << std::endl;
		#endif
	}

	template<typename _Matrix, typename _Communicator>
	struct Matrix_data {
		//Internal Matrix Type
		typedef _Matrix Matrix;
		typedef _Communicator Communicator;
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
		/*
 		* Below is a list of communicators
 		*/ 	
		_Communicator world;
		_Communicator row_comm;
		_Communicator col_comm;
		_Communicator ant_comm;
		_Communicator pan_comm;
	}; // struct Matrix_data

}//end namespace t10
#endif//T10_UTIL_H
