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

//PROJECT
#include "io.h"
#include "qr_algorithm_distributed.h"
#include "util.h"

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
	Matrix_data data;
	std::cout << env.processor_name() << " <----> " 
		  << data.world.rank() << std::endl;
	//read input
	po::variables_map vm;
	t10::process_args( argc, argv, vm);
	std::string filename( vm[ "input-file"].as< std::string>());
	t10::read_matrix( filename, data);
	t10::construct_communicators( data);
	t10::qr( data);
	std::cout << " did that error occur yet?  " << std::endl << std::flush;
}
