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
	//TODO: Delimiter should be passed into program as an option 
	template< typename Stream, typename Matrix, typename Communicator>
	bool read_csv( Stream & in, Matrix & M, const Communicator & world, char delimiter=','){
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
		//matrix is well formatted, now get the size of the blocks for each processor
		//1. calculate equally sized block sizes and remaining entries
		const std::size_t num_proc = world.size();
		const std::size_t p = std::sqrt(world.size());
		const std::size_t proc_id = world.rank();
		std::size_t block_size = number_of_rows/p;
		const std::size_t remaining = number_of_rows % p;
		//2. based on number of remaining entries and proc ID determine if block size needs to be increased by one
		//determine the starting line and starting column index
		//1. depending on the processor ID and number of left over entries
		//together with starting row and column and block size
		//1. size the matrix M this way
		//2. input entries into M
		
		//Step 2:
		//TODO: set variables like these 
		//use world.size() and world.rank()
		//const std::size_t first_row=;
		//const std::size_t last_row=;
		//const std::size_t firt_col=;
		//const std::size_t last_col=;
		//TODO: Allocate the output
		//M.resize( _____, _____);

		//Step 3: Read The File!
		//Go to Beginning...
		in.seekg(std::ios::beg);
		//... then go to the first_row line
		//for(std::size_t i =0; i < first_row-1; ++i) { in.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); }
		//for(std::size_t i =0; i < last_row-first_row; ++i){
		//	std::getline(in, line);
			//TODO: use a std::stringstream (for example) to tokenize the line
			//iterate over tokenized string and cast then insert results into M[i,j]
		//}
		return true;
	}

	template< typename Stream, typename Matrix, typename Communicator>
	bool read_mm( Stream & in, Matrix & M, const Communicator & World){
		std::cerr << "Not Yet Implemented Yet" << std::endl;
		return false;
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
			if (!read_csv( in, M, world)){
				std::cerr << "Error Reading CSV file" << std::endl;
				std::exit( -1);
			}
		} else if (file_ext == "mm") {
			if (!read_mm( in, M, world)){
				std::cerr << "Error Reading Matrix Market file" << std::endl;
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
} //end namespace t10
