cmake_minimum_required(VERSION 3.25)
project(gmpgen)

set(CMAKE_CXX_STANDARD 20)
link_libraries(gmpxx gmp)
add_executable(${PROJECT_NAME} main.cpp)
