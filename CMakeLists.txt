# to stop the whining
cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0003 NEW) # or cmake_policy(VERSION 2.6)

#fairly self explanatory
project (QR)

find_package(MPI REQUIRED)
find_package(Boost COMPONENTS mpi REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

### usually each subfolder has this
include_directories (
${MPI_INCLUDE_PATH}
${Boost_INCLUDE_PATH}
)

link_directories (${MPI_LIBRARIES})

add_executable(qrs qr_serial.cpp)
add_executable(qrp qr_distributed.cpp)

target_link_libraries(qrs ${MPI_LIBRARIES} ${Boost_LIBRARIES} boost_program_options)
target_link_libraries(qrp ${MPI_LIBRARIES} ${Boost_LIBRARIES} boost_program_options)
