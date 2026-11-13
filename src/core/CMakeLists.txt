cmake_minimum_required(VERSION 3.24)

# Preserve the current architecture default, but allow overriding it.
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 89 CACHE STRING "CUDA architectures")
endif()

project(qdyn LANGUAGES CXX CUDA)

find_package(CUDAToolkit REQUIRED)

# Retained for the existing test CMakeLists.txt.
set(CUDA_TOOLKIT_INCLUDE_DIR "${CUDAToolkit_INCLUDE_DIRS}")

option(QDYN_SPFP "Build in SPFP mixed-precision mode" OFF)
option(QDYN_BUILD_TESTS "Build QDyn test targets" ON)
option(QGPU_BUILD_PYTHON "Build the Python extension" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

file(GLOB QDYN_ALL_CUDA CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cu"
    "${CMAKE_CURRENT_SOURCE_DIR}/cuda/src/*.cu"
)

list(REMOVE_ITEM QDYN_ALL_CUDA
    "${CMAKE_CURRENT_SOURCE_DIR}/main.cu"
)

file(GLOB QDYN_ALL_CPP CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/common/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/common/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/common/src/parsers/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpu/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpu/src/*.cpp"
)

set(QDYN_CORE_SOURCES
    ${QDYN_ALL_CUDA}
    ${QDYN_ALL_CPP}
)

# One static core per precision. SPFP is a compile-time switch (precision.h),
# so each precision needs its own copy of the core.
function(qgpu_add_core target spfp)
    add_library(${target} STATIC ${QDYN_CORE_SOURCES})

    target_include_directories(${target} PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/common/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpu/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/cuda/include"
    )

    target_link_libraries(${target} PUBLIC CUDA::cudart)

    set_target_properties(${target} PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
        CUDA_RUNTIME_LIBRARY Shared
    )

    if(spfp)
        target_compile_definitions(${target} PUBLIC QDYN_SPFP)
    endif()

    target_compile_options(${target} PUBLIC
        $<$<COMPILE_LANGUAGE:CUDA>:--fmad=false>
        $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CUDA>>:-O3>
        $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CUDA>>:-G;-g>
        $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-g>
    )
endfunction()

# Core used by the CLI; precision follows -DQDYN_SPFP.
qgpu_add_core(qgpu_core ${QDYN_SPFP})

# Existing standalone executable.
add_executable(qdyn main.cu)
target_link_libraries(qdyn PRIVATE qgpu_core)

set_target_properties(qdyn PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RUNTIME_LIBRARY Shared
)

set(OUT_BIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../bin")

add_custom_target(move
    COMMAND ${CMAKE_COMMAND} -E make_directory "${OUT_BIN_DIR}"
    COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:qdyn> "${OUT_BIN_DIR}/qdyn"
    DEPENDS qdyn
)

# Optional Python extension.
if(QGPU_BUILD_PYTHON)
    set(PYBIND11_FINDPYTHON ON)
    find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
    find_package(pybind11 CONFIG REQUIRED)

    # One extension per precision: qgpu._core (double) and qgpu._core_spfp (SPFP).
    # The core matching QDYN_SPFP is reused; only the other precision is built extra.
    foreach(variant IN ITEMS double spfp)
        if(variant STREQUAL "spfp")
            set(module _core_spfp)
            set(is_spfp ON)
        else()
            set(module _core)
            set(is_spfp OFF)
        endif()

        if((is_spfp AND QDYN_SPFP) OR (NOT is_spfp AND NOT QDYN_SPFP))
            set(core qgpu_core)
        else()
            set(core qgpu_core_${variant})
            qgpu_add_core(${core} ${is_spfp})
        endif()

        pybind11_add_module(${module} MODULE NO_EXTRAS python/bindings.cpp)
        target_link_libraries(${module} PRIVATE ${core})
        target_compile_definitions(${module} PRIVATE QGPU_MODULE_NAME=${module})
    endforeach()
endif()

if(QDYN_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test/shake)
endif()