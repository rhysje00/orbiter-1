cmake_minimum_required(VERSION 3.6)
project(ORBITER)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_custom_target(makefile)
add_executable(ORBITER SRC/APPS/GRAPH/graph.out)