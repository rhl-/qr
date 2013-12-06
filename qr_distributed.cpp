#define REDIRECT_OUTPUT
//STL
#include <iostream>
#include <algorithm>
#include <string>

//BOOST UBLAS
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/banded.hpp> //for banded_adaptor

//BOOST PROGRAM OPTIONS
#include <boost/program_options.hpp>

//BOOST MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/timer.hpp>

//PROJECT
#include "io.h"
#include "util.h"
#include "communicators.h"
#include "qr_distributed.h"

namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;
namespace mpi = boost::mpi;

typedef mpi::communicator Communicator;

typedef ublas::vector< double > Vector;
typedef ublas::matrix< double> Matrix;
typedef ublas::zero_matrix< double> Zero_matrix;
typedef ublas::zero_vector< double> Zero_vector;
typedef typename ublas::diagonal_adaptor< Matrix> Diagonal_adapter;

typedef t10::Matrix_data< Matrix, Communicator> Matrix_data;

int main( int argc, char * argv[]){
	//initialize mpi
	mpi::environment env(argc, argv);
	mpi::timer t;
	Matrix_data data;
	std::cout << env.processor_name() << " <----> " 
		  << data.world.rank() << " " << t.elapsed()
		  << std::flush << std::endl;
	#ifdef REDIRECT_OUTPUT
	std::stringstream ss;
	ss << "output/out." << data.world.rank();
	std::ofstream out(ss.str().c_str());
        std::streambuf *cerrbuf = std::cerr.rdbuf(); //save old buf
        std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
        std::cerr.rdbuf(out.rdbuf()); //redirect std::cerr to out.txt!
        std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
	#endif
	std::cout << env.processor_name() << " <----> " 
		  << data.world.rank() << " " << t.elapsed()
		  <<  std::flush << std::endl;
	t.restart();
	data.world.barrier();
	std::cout << "barrier took: " << t.elapsed() << std::endl;	
	//read input
	t.restart();
	po::variables_map vm;
	t10::process_args( argc, argv, vm);
	std::string filename( vm[ "input-file"].as< std::string>());
	std::cout << "read matrix" << std::flush;
	t10::read_matrix( filename, data);
	std::cout << "... done " << t.elapsed() << std::endl 
		  << "building communicators" << std::flush;
	t.restart();
	t10::construct_communicators( data);
  	std::cout << "... done " << t.elapsed() << std::endl;	
	/*std::cout << data.world.rank() 
		  << " has a " 
		  << data.M.size1() << " x " << data.M.size2() 
		  << std::endl << std::flush;
	if (data.world.rank() == 0){ 
		std::cout << "--------------------" << std::endl << std::flush;
	}
	*/
	t10::parallel::qr( data);
	std::cerr << data.world.rank() << " is out of qr!" << std::endl;
	std::cerr << t10::print_matrix( data.M) << std::endl;
	std::cout << "-------------------------------" << std::endl << std::endl;
	#ifdef REDIRECT_OUTPUT
	std::cerr.rdbuf(cerrbuf); //reset to standard output again
	std::cout.rdbuf(coutbuf); //reset to standard output again
	out.close();
	#endif
	return 0;
}
