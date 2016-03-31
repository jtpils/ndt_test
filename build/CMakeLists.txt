cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(test)

find_package(PCL 1.2 REQUIRED)
find_package(OpenCV REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS} -g -std=c++11)

add_executable (test ../scripts/main.cpp ../scripts/ndt_test.cpp ../scripts/PSO.cpp)
target_link_libraries (test ${PCL_LIBRARIES} ${OpenCV_LIBS})
