cmake_minimum_required(VERSION 3.25)

project(rt-cpp)

set(CMAKE_CXX_STANDARD 20)

add_executable(rt-cpp src/main.cpp)

add_executable(rt-cpp-test src/test.cpp)
