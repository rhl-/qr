#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
int main () {
    using namespace boost::numeric::ublas;
    vector<double> v (5,99);
    vector_range<vector<double> > r (v, range (1, 3));
    for (unsigned i = 0; i < r.size (); ++ i)
        r (i) = i;
    for (unsigned i = 0; i < v.size(); ++i)
	std::cout << v[i] << std::endl;
    matrix<double> m (3, 3);
    for (unsigned i = 0; i < m.size1 (); ++ i) {
        matrix_row<matrix<double> > mr (m, i);
        for (unsigned j = 0; j < mr.size (); ++ j)
            mr (j) = 3 * i + j;
    }
    /*for (unsigned i = 0; i < m.size1(); ++i) 
	for(unsigned j = 0; j < m.size2(); ++j) 
		std::cout << m[i][j] << std::endl;
	*/
    std::cout << m << std::endl;
    matrix<double> m2 (3,3);
    axpy_prod(m,m,m2);
    std::cout << m2 << std::endl;
    matrix_row<matrix<double> > mr (m2,1);
    matrix<double> m21 = subrange(m2, 1,3, 0,1);
    matrix<double> m4 (2,1);
    matrix<double> m3 = subrange(m2, 1,3, 0,2);
    std::cout << m21.size1() << m21.size2() << std::endl;
    std::cout << m21 << std::endl;
    std::cout << m3<< std::endl;
    axpy_prod(m3, m21, m4);
    std::cout << m4 << std::endl;
}

