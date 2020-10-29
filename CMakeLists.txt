cmake_minimum_required(VERSION 3.10)
project(sigmod2014-pc-q3-debug)

set(CMAKE_CXX_STANDARD 11)

set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -Wextra -Wpedantic " )

link_libraries(graphblas lagraph)

find_package(OpenMP REQUIRED)
link_libraries(OpenMP::OpenMP_C)

add_executable(main main.c)
