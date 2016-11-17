cmake_minimum_required(VERSION 3.6)
project(hpcproject)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(SOURCE_FILES KMeans.c)
add_executable(hpcproject ${SOURCE_FILES})
set (CMAKE_SHARED_LINKER_FLAGS "-lm")
