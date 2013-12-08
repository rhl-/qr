include_directories(
${QR_SOURCE_DIR}/distributed
${QR_SOURCE_DIR}/serial
${QR_SOURCE_DIR}/io
${MPI_INCLUDE_PATH}
/home/rhl/boost/include
)

link_directories(
${MPI_LIBRARIES}
/home/rhl/boost/lib
)

add_executable(qrp qr_distributed.cpp)
add_executable(hessp hessenberg.cpp)
target_link_libraries(qrp ${MPI_LIBRARIES} boost_serialization 
			  boost_program_options boost_mpi)
target_link_libraries(hessp ${MPI_LIBRARIES} boost_serialization 
                          boost_program_options boost_mpi)
