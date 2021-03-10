cmake_minimum_required(VERSION 3.19) # Cmake version requirement
project(InterAtomicInteraction) # project name

set(CMAKE_CXX_STANDARD 17) # compiler flags -std=c++17
set(CMAKE_CXX_STANDARD_REQUIRED_ON)

add_subdirectory(src) # look into /src for more CMakeLists.txt instructions

enable_testing() # enable test-mode, invokable through 'ctests' command-line-tool

add_subdirectory(test) # look into /test for more CMakeLists.txt test-instructions
