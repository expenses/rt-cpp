cmake_minimum_required(VERSION 3.25)

project(rt-cpp)

find_package(OpenEXR REQUIRED)
#find_package(USD REQUIRED)

add_executable(rt-cpp src/main.cpp)

target_compile_features(rt-cpp PUBLIC cxx_std_23)
target_link_libraries(rt-cpp OpenEXR::OpenEXR)
target_compile_definitions(rt-cpp PUBLIC DBG_MACRO_NO_WARNING=1)
