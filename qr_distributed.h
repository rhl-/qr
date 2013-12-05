#ifndef QR_DISTRIBUTED_H
#define QR_DISTRIBUTED_H
#define QR_ITERATION_OUTPUT

//BOOST MPI
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

//BOOST UBLAS
#include <boost/numeric/ublas/banded.hpp> //for banded_adaptor
#include <boost/numeric/ublas/matrix.hpp> //for slices
#include <boost/numeric/ublas/matrix_proxy.hpp> //for slices
#include <boost/numeric/ublas/vector_proxy.hpp> //for slices

//PROJECT
#include "io.h"
#include "qr_serial.h"
#include "hessenberg_distributed.h"

namespace ublas = boost::numeric::ublas;
namespace t10 {
	namespace parallel {
	
	template< typename Matrix_data>
	void qr( Matrix_data & data){
			hessenberg( data);
			/*for (std::size_t i = M.size1(); i > 1; --i){
				Matrix_range R(M, Range (0, i), Range (0, i));
				qr_iteration( R);
			}*/
	}

} //end namespace parallel
} //end namespace t10
#endif //QR_DISTRIBUTED_H
