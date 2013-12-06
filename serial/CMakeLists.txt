include_directories (
${QR_SOURCE_DIR}/distributed
${QR_SOURCE_DIR}/serial
${QR_SOURCE_DIR}/io
/home/rhl/boost/include
)

link_directories (
/home/rhl/boost/lib
)

add_executable(qrs qr_serial.cpp)
target_link_libraries(qrs ${Boost_LIBRARIES} boost_program_options)
