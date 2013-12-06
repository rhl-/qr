# to stop the whining
cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0003 NEW) # or cmake_policy(VERSION 2.6)

#fairly self explanatory
project (QR)

find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_subdirectory(serial)
add_subdirectory(distributed)
add_subdirectory(io)
