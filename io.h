#ifndef T10_IO_H
#define T10_IO_H
#include <boost/numeric/ublas/io.hpp> //io
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include "mmio.h"

//BOOST MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

//BOOST PROGRAM OPTIONS
#include <boost/program_options.hpp>

namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;
namespace t10 {
	
        std::pair<std::size_t, std::size_t> id_to_index(const std::size_t & proc_id, const std::size_t & p) {
                return std::make_pair(proc_id / p, proc_id % p);
        }

        std::size_t index_to_id(const std::size_t & block_row, const std::size_t & block_col, const std::size_t & p) {
                return block_row * p + block_col;
        }
	//TODO: Delimiter should be passed into program as an option 
	template< typename Stream, typename Matrix_data, typename Communicator>
	bool read_csv( Stream & in, Matrix_data & data, const Communicator & world, char delimiter=','){
		typedef typename Matrix_data::Matrix Matrix;
		typedef typename std::pair<std::size_t, std::size_t> Block;

		Matrix & M = data.M;
		std::string line;

		//Step 0: Get first line
		std::getline(in, line);
		std::size_t number_of_rows=1;
		const std::size_t number_of_columns=std::count(line.begin(), line.end(), delimiter)+1;

		//Step 1: Determine the number of lines in the file
		//Also verify that it is well formatted.
		while (std::getline(in, line)){
			if(number_of_columns!=std::count(line.begin(), line.end(), delimiter)+1){
				std::cerr << "CSV File is invalid" << std::endl;
				return false;
			}
			++number_of_rows;
		}
		if (number_of_rows != number_of_columns){
			std::cerr << "Matrix is not square" << std::endl;
			return false;
		}
		data.n = number_of_rows;
		//matrix is well formatted, now get the size of the blocks for each processor
		//1. calculate equally sized block sizes and remaining entries
		const std::size_t num_proc = world.size();
		// for now ignore the fact that number of processors is not perfect square
		const std::size_t p = std::sqrt(world.size());
		const std::size_t proc_id = world.rank();
		const std::size_t avg_block_size = number_of_rows / p;
		// row size of block before resizing
		std::size_t block_size1 = avg_block_size;
		// col size of block before resizing
		std::size_t block_size2 = avg_block_size;
		// number of left-over elements
		const std::size_t remaining = number_of_rows % p;
		// row and column index of the block for this processor, from 0 to p-1, this is invariant
		Block block = id_to_index(proc_id, p);
		
		const std::size_t block_row = block.first;
		const std::size_t block_col = block.second;
		data.block_row = block_row;
		data.block_col = block_col;
		// if remaining was > 0, this block size increases by 1 if the block is within the first remaining x remaining 
		// blocks, as we distribute remaining entries among the first submatrix of blocks

		block_size1 += block_row < remaining;
		block_size2 += block_col < remaining;
		// number of lines to skip to get to this block, sum of the block sizes of the blocks "above" it
		// add the lesser of remaining number of elements and block_row index, because that many remaining elements
		// have been distributed to the blocks "above" it
		const std::size_t first_row = block_row * avg_block_size + std::min(remaining, block_row);
		// similarly for columns
		const std::size_t first_col = block_col * avg_block_size + std::min(remaining, block_col);
		data.first_row = first_row;
		data.first_col = first_col;
		data.last_row = first_row + block_size1;
		data.last_col = first_col + block_size2;
		// now set the matrix size
		M.resize(block_size1, block_size2);

		in.clear();
		in.seekg(0, std::ios::beg);
		//Step 3: Read The File!
		//Go to Beginning...
		//... then go to the first_row line
		for(std::size_t i =0; i < first_row; ++i) { std::getline(in, line); }
		for(std::size_t i =0; i < block_size1; ++i){
			std::getline(in, line);
			for (std::size_t j = 0; j < first_col; ++j){ 
				const std::size_t found = line.find_first_of(delimeter);
				line = line.substr( found+1);
			}

			for( std::size_t j = 0; j < block_size2; ++j){
				const std::size_t found = line.find_first_of(delimeter);
				M(i,j) = atof( line.substr(0, found).c_str());
				line = line.substr( found+1);
			}
		}
		return true;
	}

	template< typename Stream, typename Matrix, typename Communicator>
	bool read_mm( Stream & in, Matrix & M, const Communicator & World){
		std::cerr << "Not Yet Implemented Yet" << std::endl;
		return false;
	}
	template< typename Stream, typename Matrix, typename Communicator>
	bool read_mat( Stream & in, Matrix & M, const Communicator & World){
		std::cerr << "Not Yet Implemented Yet" << std::endl;
		return false;
	}
	template< typename String, typename Matrix, typename Communicator>
	void read_matrix( const String & filename, Matrix & M, 
			  const Communicator & world){
		const std::string file_ext(filename.substr(filename.find_last_of(".") + 1));
		std::ifstream in(filename.c_str());
		if (!in.good()){ 
			std::cerr << "Error Opening " << filename  << std::endl;
			std::exit( -1);
		}
		//TODO: Error checking if file doesn't open.
		if (file_ext == "csv") {
			if (!read_csv( in, M, world)){
				std::cerr << "Error Reading CSV file" << std::endl;
				std::exit( -1);
			}
		} else if (file_ext == "mm") {
			if (!read_mm( in, M, world)){
				std::cerr << "Error Reading Matrix Market file" << std::endl;
				std::exit( -1);
			}
		} else if (file_ext == "mat") {
			if (!read_mat( in, M, world)){
				std::cerr << "Error Reading .Mat file" << std::endl;
				std::exit( -1);
			}
		} else {
			std::cerr << "File Format Not Supported" << std::endl;
			std::exit( -1);
		}
	}
	template<typename Variable_map>
	void process_args( int & argc, char *argv[],Variable_map & vm){
	  //parse command line options
	  po::options_description desc( "Usage: tp [options] input-file");
	  desc.add_options()
	  ( "help", "Display this message")
	  ( "input-file", "input matrix file to parse");
	  po::positional_options_description p;
	  p.add( "input-file",1);
	  po::store( po::command_line_parser( argc,argv)
	            .options( desc)
	            .positional( p)
	            .run(),
	             vm);
	  po::notify( vm);
	  if ( vm.count( "help")){
	        std::cout << desc << std::endl;
	        std::exit( 0);
	  }
	  if ( vm.count( "input-file") != 1){
	        std::cout << desc << std::endl;
	        std::exit( -1);
	  }
	}

	template< typename Matrix>
	std::string print_matrix( const Matrix & M, std::size_t p=7){
		std::string str;
		std::stringstream ss;
		ss.setf( std::ios::fixed, std:: ios::floatfield );
		ss << std::setw( p) << std::setprecision( p) <<  M;
		str = ss.str();
		//kill size thing
		std::size_t lefts = str.find(std::string("["));
		std::size_t rights = str.find(std::string("]"));
		str.replace(lefts,rights-lefts+1,"");
		
		//start with new line	
		std::size_t found = str.find(std::string("(("));
		str.replace(found,2,std::string("\n"));
		//end with new line
		found = str.find(std::string("))"));
		str.replace(found,2,std::string("\n"));
	
		//end of each row
		found = str.find(std::string("),"));
		while( found != std::string::npos){
			str.replace(found,2,std::string("\n"));
			found = str.find(std::string("),"));
		}
		//remove beginning of row
		found = str.find(std::string("("));
		while( found != std::string::npos){
			str.replace(found,1,std::string(""));
			found = str.find(std::string("("));
		}
		//replace , with tab
		found = str.find(std::string(","));
		while( found != std::string::npos){
			str.replace(found,1,std::string("\t"));
			found = str.find(std::string(","));
		}
		
		return str;
	}



	
	template< typename Matrix>
	std::string print_matrix( const Matrix & M, bool matlab){
		std::string str;
		std::stringstream ss;
		ss << M;
		str = ss.str();
		std::size_t lefts = str.find(std::string("["));
		std::size_t rights = str.find(std::string("]"));
		str.replace(lefts,rights-lefts+1,"");
		std::size_t found = str.find(std::string("(("));
		str.replace(found,2,std::string("["));
		found = str.find(std::string("))"));
		str.replace(found,2,std::string("]"));
		found = str.find(std::string("),"));
		while( found != std::string::npos){
			str.replace(found,2,std::string(";"));
			found = str.find(std::string("),"));
		}

		found = str.find(std::string("("));
		while( found != std::string::npos){
			str.replace(found,1,std::string(""));
			found = str.find(std::string("("));
		}


		return str;	
	}
}
#endif
 //end namespace t10
