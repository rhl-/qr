#ifndef T10_UTIL_H                                          
#define T10_UTIL_H                                          
//#define UTIL_DEBUG 
//Boost MPI
#include <boost/mpi.hpp>
#include <algorithm>
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
		typedef typename std::vector< Communicator> Vector_comm;	
		const Communicator & world = data.world;

		const std::size_t id = world.rank();
		const std::size_t world_size = world.size();
		const std::size_t row_length = std::sqrt(world_size);
		const std::size_t r = row_id( id, row_length);
		const std::size_t c = col_id( id, row_length);
		const std::size_t diag1 = std::min(r,c);
		const std::size_t diag2 = std::max(r,c);
		const std::size_t anti_m = std::min(r+c, row_length-1);
		const std::size_t anti_size = 2*anti_m-(r+c) + 1;
		
		data.partner = index_to_id( c,r,  row_length);

		Vector & ant = data.ant;
		Vector & row = data.row;
		Vector & col = data.col;
		Vector & pan = data.pan;
		Vector & l_col = data.l_col;
		Vector & r_row = data.row;
		Vector & p_row = data.p_row;
		Vector & p_col = data.p_col;
		Vector & s_col = data.s_col;
		
		//communicators involved in right householder multiplication (nvolves partner col and own row)
		Vector_comm & right_comm_vec = data.right_comm_vec;
		//comm involved in left mult (needs own col and partner row)
		Vector_comm & left_comm_vec = data.left_comm_vec;
		//communicators along row of this processor (for propagating right multiplication)
		Vector_comm & row_comm_vec = data.row_comm_vec;
		//communicators along col of this processor (for propagating left multiplication)
		Vector_comm & col_comm_vec = data.col_comm_vec;

		ant.resize( anti_size, id);
		row.resize( row_length, id);
		col.resize( row_length, id);
		pan.resize( 2*(row_length-diag1)-1, id);
		l_col.resize( row_length - diag1, id);
		s_col.resize( l_col.size()-1, id);
		r_row.resize( row_length - diag1, id);
		p_row.resize( row_length, data.partner);
		p_col.resize( row_length, data.partner);
		
		//vector of all panels of communicators
		Vector_comm all_panels(row_length);		
		//loop over diag entries from 0 to p-1
		diag = 0;
		for (Iterator k = all_panels.begin(); k != all_panels.end(); ++k) {
			Vector pan_j(2*(row_length-diag)-1, id);
			Iterator end = pan_j.begin()+(row_length-diag);
                	std::size_t a = diag;
                	std::size_t b = diag;
                	for (Iterator i = pan_j.begin(); i != end; ++i,++a){
                        	*i = index_to_id(a,diag, row_length);
                	}

                	for(Iterator i = end; i != pan_j.end(); ++i){
                	        *i = index_to_id(diag,++b, row_length);
        	        }
			mpi::group pan_j_group = world.group().include( pan_j.begin(),
                                                                pan_j.end());
			++diag;
		}
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

		//create lower col indices
		for (Iterator j = l_col.begin(); j != l_col.end(); ++j) {
			const std::size_t dst = std::distance(l_col.begin(), j);
			const std::size_t k = dst + diag1;
			*j = k*row_length + c;
		}

		//create right row indices
		for (Iterator i = r_row.begin(); i != r_row.end(); ++i) {
			const std::size_t dst = std::distance(r_row.begin(), i);
			const std::size_t k = dst + diag1;
			*i = r*row_length + k;
		}

		//create partner proc col and row by observing that the only 
		//difference is swapping c and r
		j = p_col.begin();
                for(Iterator i = p_row.begin(); i != p_row.end(); ++i, ++j){
                        const std::size_t k = std::distance(p_row.begin(), i);
                        *i = c*row_length + k;
                        *j = k*row_length + r;
                }

		std::sort (row.begin(),   row.end()); 
		std::sort (col.begin(),   col.end()); 
		std::sort (ant.begin(),   ant.end()); 
		std::sort (pan.begin(),   pan.end()); 
		std::sort (l_col.begin(), l_col.end()); 
		std::sort (r_row.begin(), r_row.end()); 
		std::sort (p_row.begin(), p_row.end()); 
		std::sort (p_col.begin(), p_col.end());
		std::copy (l_col.begin()+1,l_col.end(),s_col.begin());
		#ifdef UTIL_DEBUG
			std::cout << "Processor: " << world.rank() <<  std::endl
				  << "Partner: " << data.partner << std::endl
				  << "row: " << row << std::endl 
				  << "col: " << col << std::endl 
				  << "ant: " << ant << std::endl
				  << "pan: " << pan << std::endl 
                                  << "Lower col: " << l_col << std::endl
				  << "Right row: " << r_row << std::endl
				  << "Partner row: " << p_row << std::endl
				  << "Partner col: " << p_col << std::endl;
		#endif
		mpi::group row_group = world.group().include( row.begin(), 
								row.end());
		mpi::group col_group = world.group().include( col.begin(), 
								col.end());
		mpi::group ant_group = world.group().include( ant.begin(), 
								ant.end());
		mpi::group pan_group = world.group().include( pan.begin(), 
								pan.end());
		mpi::group l_col_group = world.group().include( l_col.begin(),
                                                                l_col.end());
		mpi::group s_col_group = world.group().include( s_col.begin(),
                                                                s_col.end());
		mpi::group r_row_group = world.group().include( r_row.begin(),
                                                                r_row.end());
		mpi::group p_col_group = world.group().include( p_col.begin(),
                                                                p_col.end());
		mpi::group p_row_group = world.group().include( p_row.begin(),
                                                                p_row.end());
		data.row_comm = mpi::communicator(world, row_group);
		data.col_comm = mpi::communicator(world, col_group);
		data.ant_comm = mpi::communicator(world, ant_group);
		data.pan_comm = mpi::communicator(world, pan_group);
		data.r_row_comm = mpi::communicator(world, r_row_group);
                data.l_col_comm = mpi::communicator(world, l_col_group);
                data.s_col_comm = mpi::communicator(world, s_col_group);
                data.p_col_comm = mpi::communicator(world, p_col_group);
                data.p_row_comm = mpi::communicator(world, p_row_group);

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
		
                          << data.l_col_comm.rank() << " of "
                          << data.l_col_comm.size() << " in l_col_comm"
                          << std::endl

                          << data.r_row_comm.rank() << " of "
                          << data.r_row_comm.size() << " in r_row_comm"
                          << std::endl

                          << data.p_col_comm.rank() << " of "
                          << data.p_col_comm.size() << " in p_col_comm"
                          << std::endl

                          << data.p_row_comm.rank() << " of "
                          << data.p_row_comm.size() << " in p_row_comm"
                          << std::endl
			  << "-------------------------" 
			  << std::endl;
		#endif
	}

	template<typename _Matrix, typename _Communicator>
	struct Matrix_data {
		typedef _Matrix Matrix;
		typedef _Communicator Communicator;
		typedef typename std::vector< std::size_t> Vector;
		typedef typename std::vector< _Communicator> Vector_comm;

		bool above() const { return block_row <= block_col;}
		bool below() const { return block_col <= block_row;}
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
		Vector ant;
		Vector row;
		Vector col;
		Vector pan;
		// processor below
		Vector l_col;
		//sub column
		Vector s_col;
		// proc to the right of
		Vector r_row;
		// proc on partner proc row
		Vector p_row;
		// proc on partner proc col
		Vector p_col;
		/*
 		* Below is a list of communicators
 		*/ 	
		_Communicator world;
		_Communicator row_comm;
		_Communicator col_comm;
		_Communicator ant_comm;
		_Communicator pan_comm;
		_Communicator r_row_comm;
                _Communicator l_col_comm;
                _Communicator s_col_comm;
                _Communicator p_row_comm;
                _Communicator p_col_comm;
		Vector_comm right_comm_vec;
		Vector_comm left_comm_vec;
		Vector_comm col_comm_vec;
		Vector_comm row_comm_vec;
	}; // struct Matrix_data

}//end namespace t10
#endif//T10_UTIL_H
