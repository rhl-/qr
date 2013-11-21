#ifndef QR_ALGORITHM_H                                          
#define QR_ALGORITHM_H                                          
#define QR_ITERATION_OUTPUT  
namespace t10 {
template<typename _Matrix>
        struct Matrix_data {
        typedef typename _Matrix Matrix;
        Matrix M;
        std::size_t n;
        std::size_t first_row;
        std::size_t last_row;
        std::size_t first_col;
        std::size_t last_col;
        std::size_t block_row;
        std::size_t block_col;
}; // struct Matrix_data
}//end namespace t10
#endif//QR_ALGORITHM_H
