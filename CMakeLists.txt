cmake_minimum_required(VERSION 3.25)
project(gmpgen)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP)

link_libraries(gmpxx gmp)
link_libraries(OpenMP::OpenMP_CXX)

add_executable(${PROJECT_NAME} main.cpp gmpgen.cpp)

