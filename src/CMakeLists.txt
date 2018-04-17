# At LEAST 2.8 but newer is better
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(ecare)

cmake_minimum_required(VERSION 2.8.12)
project(ecare)

add_subdirectory(../modules/dlib ../build/dlib_build)
include_directories(../include)
add_executable(ecare ecare.cpp ellipse.cpp)
target_link_libraries(ecare dlib::dlib)