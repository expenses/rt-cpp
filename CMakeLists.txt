cmake_minimum_required(VERSION 3.25)

project(rt-cpp)

find_package(OpenEXR REQUIRED)

include(external/usd/pxrConfig.cmake)

file(GLOB SRC_FILES src/*.cpp)
add_executable(rt-cpp ${SRC_FILES})

include_directories(
    external
    external/usd/include
)
target_compile_features(rt-cpp PUBLIC cxx_std_20)
target_link_libraries(rt-cpp OpenEXR::OpenEXR ${PXR_LIBRARIES})
target_compile_definitions(rt-cpp PUBLIC DBG_MACRO_NO_WARNING=1)

# Only for select headers
target_precompile_headers(
    rt-cpp PRIVATE
    src/external.h src/util.h src/sampling.h
)
