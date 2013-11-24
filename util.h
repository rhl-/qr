#ifndef T10_UTIL_H                                          
#define T10_UTIL_H                                          
//#define UTIL_DEBUG 
//Boost MPI
#include <boost/mpi.hpp>
#include <algorithm>
namespace mpi = boost::mpi;
namespace t10 {

	template< typename Communicator, typename Vector> 
	boost::mpi::group create_group( const Communicator & world, 
						   const Vector & v){ 
		return world.group().include(v.begin(), v.end()); 
	}
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
	//TODO: Non-exported functionality
	template< typename Vector>
	void create_panel_vector( Vector & panel, const std::size_t k, 
						  const std::size_t row_length){
		typedef typename Vector::iterator Iterator;
		Iterator end = panel.begin()+(row_length-k);
        	std::size_t a = k, b = k;
        	for (Iterator i = panel.begin(); 
			      i != end; ++i,++a){
        		*i = index_to_id( a, k, row_length);
        	}

        	for(Iterator i = end; i != panel.end(); ++i){
        	        *i = index_to_id( k, ++b, row_length);
        	}
	}

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
		const std::size_t diag1 = std::min(r,c);
		const std::size_t diag2 = std::max(r,c);
		const std::size_t anti_m = std::min(r+c, row_length-1);
		const std::size_t anti_size = 2*anti_m-(r+c) + 1;
		
		data.partner = index_to_id( c,r,  row_length);


		Vector row( row_length, id);
		Vector col( row_length, id);
		Vector l_col( row_length - diag1, id);
		Vector s_col( l_col.size()-1, id);
		Vector r_row( row_length - diag1, id);
		Vector p_row( row_length, data.partner);
		Vector p_col( row_length, data.partner);

		//create row and column mpi group indices
		//for proc and its partner
		Iterator j = col.begin();
		Iterator k = p_col.begin();
		Iterator l = p_row.begin();
		for(Iterator i = row.begin(); i != row.end(); ++i, ++j,++k,++l){
			const std::size_t idx = std::distance(row.begin(), i);
			*i = r*row_length + idx;
			*k = c*row_length + idx;
			*j = idx*row_length + c;
                        *l = idx*row_length + r;
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

		std::sort (row.begin(),   row.end()); 
		std::sort (col.begin(),   col.end()); 
		std::sort (l_col.begin(), l_col.end()); 
		std::sort (r_row.begin(), r_row.end()); 
		std::sort (p_row.begin(), p_row.end()); 
		std::sort (p_col.begin(), p_col.end());
		std::copy (l_col.begin()+1,l_col.end(),s_col.begin());
	
		mpi::group row_group   = t10::create_group (world, row);
		mpi::group col_group   = t10::create_group (world, col);
		mpi::group l_col_group = t10::create_group (world, l_col);
		mpi::group s_col_group = t10::create_group (world, s_col);
		mpi::group r_row_group = t10::create_group (world, r_row);
		mpi::group p_col_group = t10::create_group (world, p_col);
		mpi::group p_row_group = t10::create_group (world, p_row);

		data.r_row_comm = mpi::communicator(world, r_row_group);
                data.l_col_comm = mpi::communicator(world, l_col_group);
                data.s_col_comm = mpi::communicator(world, s_col_group);
                data.p_col_comm = mpi::communicator(world, p_col_group);
                data.p_row_comm = mpi::communicator(world, p_row_group);

		//communicators involved in right householder multiplication 
		//(nvolves partner col and own row)
		Vector_comm  & right_comm = data.right_comm;
		//comm involved in left mult (needs own col and partner row)
		Vector_comm  & left_comm = data.left_comm;
		//communicators along row of this processor 
		//(for propagating right multiplication)
		Vector_comm  & row_comm = data.row_comm;
		//communicators along col of this processor 
		//(for propagating left multiplication)
		Vector_comm  & col_comm = data.col_comm;

		//first union these groups
		mpi::group right_group = p_col_group | row_group;
		mpi::group left_group = p_row_group | col_group;
		mpi::group row_comm_group = row_group;
		mpi::group col_comm_group = col_group;

		//put them in the appropriate places
		right_comm.push_back( mpi::communicator( world, right_group));
		left_comm.push_back( mpi::communicator( world, left_group));
		row_comm.push_back( mpi::communicator( world, row_comm_group));
		left_comm.push_back( mpi::communicator( world, col_comm_group));
		
		//now moving down the "diagonal" remove stuff from each
		//group making communicators	
		Vector indices;
		indices.reserve( row_length);
		for (std::size_t k = 0; k < std::min(row_length-1,id); ++k) {
			indices.push_back( k);
			row_comm_group.exclude( indices.begin(), indices.end());
			col_comm_group.exclude( indices.begin(), indices.end());
			row_comm.push_back( mpi::communicator( world, 
							       row_comm_group));
			left_comm.push_back( mpi::communicator( world, 
							       col_comm_group));
			
			Vector panel(2*(row_length-k)-1);
			t10::create_panel_vector( panel, k, row_length);
			mpi::group panel_group = world.group().include( 
								panel.begin(),
                                                                panel.end());
			//set difference
			left_group = left_group - panel_group;
			right_group = right_group - panel_group;
			right_comm.push_back( mpi::communicator( world, 
								 right_group));
			left_comm.push_back( mpi::communicator( world, 
								left_group));
		}
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
		
		/*
 		* Below is a list of communicators
 		*/ 	
		_Communicator world;
		_Communicator ant_comm;
		_Communicator r_row_comm;
                _Communicator l_col_comm;
                _Communicator s_col_comm;
                _Communicator p_row_comm;
                _Communicator p_col_comm;
		Vector_comm right_comm;
		Vector_comm left_comm;
		Vector_comm col_comm;
		Vector_comm row_comm;
	}; // struct Matrix_data

}//end namespace t10
#endif//T10_UTIL_H
