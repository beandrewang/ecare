# At LEAST 2.8 but newer is better

cmake_minimum_required(VERSION 2.8.12)
project(ecare)

add_subdirectory(../modules/armadillo-code ../build/arma)
include_directories(../include)
include_directories(../modules/armadillo-code/include)
add_executable(ecare ecare.cpp ellipse.cpp)
add_executable(matrix_test test/matrix_test.cpp)
target_link_libraries(ecare armadillo)
target_link_libraries(matrix_test armadillo)
set_property(TARGET ecare PROPERTY CXX_STANDARD 14)
set_property(TARGET matrix_test PROPERTY CXX_STANDARD 14)
