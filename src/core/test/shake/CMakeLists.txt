add_executable(shake_cpu_gpu_compare EXCLUDE_FROM_ALL
  shake_cpu_gpu_compare.cu
  ${QDYN_CORE_SOURCES}
)

target_include_directories(shake_cpu_gpu_compare
  PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${CUDA_TOOLKIT_INCLUDE_DIR}"
    "${PROJECT_SOURCE_DIR}/cuda/include"
    "${PROJECT_SOURCE_DIR}/common/include"
    "${PROJECT_SOURCE_DIR}/cpu/include"
)

if(QDYN_SPFP)
  target_compile_definitions(shake_cpu_gpu_compare PRIVATE QDYN_SPFP)
endif()

target_compile_options(shake_cpu_gpu_compare
  PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--fmad=false>
    $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CUDA>>:-O3>
    $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CUDA>>:-G;-g>
    $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-g>
)

set_target_properties(shake_cpu_gpu_compare PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_test(NAME shake_cpu_gpu_compare COMMAND shake_cpu_gpu_compare)
set_tests_properties(shake_cpu_gpu_compare PROPERTIES SKIP_RETURN_CODE 77)

add_executable(cuda_serial_shake_regression
  cuda_serial_shake_regression.cu
  ${PROJECT_SOURCE_DIR}/cuda/src/cuda_serial_shake.cu
)

target_include_directories(cuda_serial_shake_regression
  PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${CUDA_TOOLKIT_INCLUDE_DIR}"
    "${PROJECT_SOURCE_DIR}/cuda/include"
    "${PROJECT_SOURCE_DIR}/common/include"
    "${PROJECT_SOURCE_DIR}/cpu/include"
)

target_compile_options(cuda_serial_shake_regression
  PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--fmad=false>
    $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CUDA>>:-O3>
    $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CUDA>>:-G;-g>
)

set_target_properties(cuda_serial_shake_regression PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_test(NAME cuda_serial_shake_regression COMMAND cuda_serial_shake_regression)
set_tests_properties(cuda_serial_shake_regression PROPERTIES SKIP_RETURN_CODE 77)
