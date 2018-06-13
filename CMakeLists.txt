cmake_minimum_required(VERSION 2.8)
project(KDTree CXX)

set(CMAKE_CXX_COMPILER             "/usr/bin/g++")
set(CMAKE_CXX_FLAGS                "-Wall -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG          "-g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

set(CMAKE_AR      "/usr/bin/gcc-ar")
set(CMAKE_LINKER  "/usr/bin/ld")
set(CMAKE_NM      "/usr/bin/gcc-nm")
set(CMAKE_OBJDUMP "/usr/bin/objdump")
set(CMAKE_RANLIB  "/usr/bin/gcc-ranlib")

set(CMAKE_BUILD_TYPE Debug)

include_directories(.)
add_library(KDTree SHARED "KDTree.cpp")

#add_executable(test_kdtree "test_kdtree.cpp")
#target_link_libraries(test_kdtree KDTree)
