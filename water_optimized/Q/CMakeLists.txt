set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
cmake_minimum_required(VERSION 3.18)
project(MyCudaProject LANGUAGES CXX CUDA)

set(CMAKE_CUDA_ARCHITECTURES 86)

enable_language(CUDA)

file(GLOB_RECURSE SRC_FILES
    ${CMAKE_SOURCE_DIR}/src/core/*.cu
    ${CMAKE_SOURCE_DIR}/src/core/*.h
)

add_executable(MyCudaApp ${SRC_FILES})

set_target_properties(MyCudaApp PROPERTIES
    CXX_STANDARD 17
    CUDA_STANDARD 17
)

target_compile_options(MyCudaApp PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
)

target_include_directories(MyCudaApp PRIVATE ${CMAKE_SOURCE_DIR}/src/core)
