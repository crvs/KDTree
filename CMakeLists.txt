cmake_minimum_required(VERSION 2.8)
project(KDTree CXX)

set(CMAKE_CXX_COMPILER             "/usr/bin/g++")
set(CMAKE_CXX_FLAGS                "-Wall -std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG          "-g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

set(CMAKE_AR      "/usr/bin/gcc-ar")
set(CMAKE_LINKER  "/usr/bin/ld")
set(CMAKE_NM      "/usr/bin/gcc-nm")
set(CMAKE_OBJDUMP "/usr/bin/objdump")
set(CMAKE_RANLIB  "/usr/bin/gcc-ranlib")

set(CMAKE_BUILD_TYPE RELEASE)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

include_directories(.)
add_library(KDTree SHARED KDTree.cpp)
target_include_directories(KDTree PUBLIC .)
add_subdirectory(tests)
