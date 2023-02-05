cmake_minimum_required(VERSION 3.25)

project(rt-cpp)

find_package(OpenEXR REQUIRED)
find_package(embree REQUIRED)
find_package(OpenImageIO REQUIRED)

include(external/usd2/pxrConfig.cmake)

file(GLOB SRC_FILES src/*.cpp)
add_executable(rt-cpp ${SRC_FILES})

include_directories(
    external
    external/usd2/include
)
target_compile_features(rt-cpp PUBLIC cxx_std_20)
target_link_libraries(rt-cpp OpenEXR::OpenEXR ${PXR_LIBRARIES} embree OpenImageIO::OpenImageIO)
target_compile_definitions(rt-cpp PUBLIC DBG_MACRO_NO_WARNING=1)

# Only for select headers
target_precompile_headers(
    rt-cpp PRIVATE
    src/external.hpp src/util.hpp src/sampling.hpp
)

target_compile_options(rt-cpp PRIVATE
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
        -Weverything -Wno-c++98-compat -Wno-pre-c++20-compat-pedantic -Wno-c++98-compat-pedantic -Wno-shadow-field-in-constructor -Wno-shadow -Wno-missing-prototypes -Wno-header-hygiene -Wno-padded
    >
)
