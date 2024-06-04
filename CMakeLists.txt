cmake_minimum_required(VERSION 3.6)
project(KDTree CXX)

#set(CMAKE_CXX_COMPILER             "/usr/bin/g++")
set(CMAKE_CXX_FLAGS                "-Wall -Wextra -Werror -std=c++17")
set(CMAKE_CXX_FLAGS_DEBUG          "-g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
	message(WARNING "Build type unset using default")
endif()
message(STATUS "Using CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

set(CMAKE_AR      "/usr/bin/gcc-ar")
set(CMAKE_LINKER  "/usr/bin/ld")
set(CMAKE_NM      "/usr/bin/gcc-nm")
set(CMAKE_OBJDUMP "/usr/bin/objdump")
set(CMAKE_RANLIB  "/usr/bin/gcc-ranlib")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

include_directories(.)
add_library(KDTree SHARED KDTree.cpp)
target_include_directories(KDTree PUBLIC .)

include(CTest)
enable_testing()
add_subdirectory(tests)

file(GLOB ALL_SOURCE_FILES *.cpp *.hpp **/*.cpp **/*.hpp)
add_custom_target(check_clang_tidy clang-tidy-14 -p ${CMAKE_CURRENT_BINARY_DIR} --config-file=${CMAKE_CURRENT_SOURCE_DIR}/.clang-tidy.yml ${ALL_SOURCE_FILES})
add_custom_target(check_clang_format clang-format-14 --style=file:${CMAKE_CURRENT_SOURCE_DIR}/.clang-format.yml --dry-run -Werror ${ALL_SOURCE_FILES})
add_custom_target(run_clang_format clang-format-14 --style=file:${CMAKE_CURRENT_SOURCE_DIR}/.clang-format.yml -i ${ALL_SOURCE_FILES})
