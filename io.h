#include <boost/numeric/ublas/io.hpp> //io
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>

//BOOST MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

//BOOST PROGRAM OPTIONS
#include <boost/program_options.hpp>

namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;
namespace t10 {
	template< typename Stream, typename Matrix, typename Communicator>
	bool read_csv( Stream & in, Matrix & M, const Communicator & World){
		const std::size_t proc_id = world.rank();
		const std::size_t num_proc = world.size();
		int N = 0;
		//Make sure to go top down.
		while(std::getline(ifs,line,',')) {
			std::stringstream lineStream(line);
			std::string token;
			while(lineStream >> token) { N++; }
		}
			int n = std::sqrt(N);
			if (N > n*n) {
				std::cout << "File has a non-square matrix" << std::endl;
				return;
			}
		return true;
	}
	template< typename Stream, typename Matrix, typename Communicator>
	bool read_mm( Stream & in, Matrix & M, const Communicator & World){
		const std::size_t proc_id = world.rank();
		const std::size_t num_proc = world.size();
		return true;
	}
	
	template< typename String, typename Matrix, typename Communicator>
	void read_matrix( const String & filename, Matrix & M, const Communicator & world){
		//0. open file (std::istream)
		//1. determine file type
			//a) determined by extension
			//b) valid extensions
				//csv
				//matrixmarket
		//2. determine matrix size k
			//take N # of rows in original matrix, divide by sqrt( num_processors)
		//3. allocate output memory i.e. k x k matrix
				//M.resize(k,k)
		//4. read appropriate part of input into output
				//seek to place in file and read elt by elt

		const std::string file_ext(filename.substr(filename.find_last_of(".") + 1));
		std::ifstream in(filename.c_str());
		if (!in.good()){ 
			std::cerr << "Error Opening " << filename  << std::endl;
			std::exit( -1);
		}
		//TODO: Error checking if file doesn't open.
		if (file_ext == "csv") {
			if (!read_csv( ifs, M, world)){
				std::cerr << "Error Reading CSV file" << std::endl;
				std::exit( -1);
			}
		} else if (file_ext == "mm") {
			if (!read_mm( ifs, M, world)){
				std::cerr << "Error Reading Matrix Market file" << std::endl;
				std::exit( -1);
			}
		} else {
			std::cerr << "File Format Not Supported" << std::endl;
			std::exit( -1);
		}
		/*
		//Get file extension
		std::string line;
		if (fileext == "csv") {
			int N = 0;
			 while(std::getline(ifs,line,',')) {
				std::stringstream lineStream(line);
				std::string token;
				while(lineStream >> token) {
					N++;
				}
			}
			int n = std::sqrt(N);
			if (N > n*n) {
				std::cout << "File has a non-square matrix" << std::endl;
				return;
			}
		}	
		else if (fileext == "mm") return;
			// matrix market format
		
		//out of world.size() total processors	
		//initially M is not the appropriate size
		// 1. open the file in filename
		// 2. we need to be able to know the matrix size in advance, say it is n x n
		// 3. M.resize(block_size(n),block_size(n))  
		// 4. read in data into M
		*/
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
} //end namespace t10
