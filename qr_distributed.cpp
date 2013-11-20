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
#include "qr_algorithm_serial.h"

namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;
namespace mpi = boost::mpi;

typedef ublas::vector< double > Vector;
typedef ublas::matrix< double> Matrix;
typedef ublas::zero_matrix< double> Zero_matrix;
typedef ublas::zero_vector< double> Zero_vector;
typedef typename ublas::diagonal_adaptor< Matrix> Diagonal_adapter;

int main( int argc, char * argv[]){
	//initialize mpi
	mpi::environment env(argc, argv);
  	mpi::communicator world;
	//read input
	po::variables_map vm;
	t10::process_args( argc, argv, vm);
	std::string filename( vm[ "input-file"].as< std::string>());
	Matrix M;
	t10::read_matrix( filename, M, world);
	std::cout << t10::print_matrix( M) << std::endl;
//	t10::qr(M, world);
}
