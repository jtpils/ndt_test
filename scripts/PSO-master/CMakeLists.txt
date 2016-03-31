cmake_minimum_required(VERSION 2.8)
ADD_DEFINITIONS(
  -O2
  -Wall
  -g
  -std=c++11
)
project( PSO )

add_executable( PSO main.cpp PSO.cpp )
target_link_libraries( PSO )