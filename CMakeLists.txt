#include_directories(${PROJECT_SOURCE_DIR}/include)
#link_directories(${PROJECT_SOURCE_DIR}/lib)
cmake_minimum_required(VERSION 2.6)
project(Voxelizer)

IF( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -DNDEBUG -mavx -Wmaybe-uninitialized ")
    set(CMAKE_EXE_LINKER_FLAGS "-DNDEBUG -mavx")
ELSEIF ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -DNDEBUG")
ENDIF()

file(GLOB SRC_LIST "*.cpp")
file(GLOB HEADER_LIST "*.h")

file(GLOB TEST_SRC "./test/*.cpp")

add_executable(voxelizer ${SRC_LIST} ${HEADER_LIST})
add_executable(test ${TEST_SRC})
