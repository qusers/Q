# Host-only unit test of the Gapsys soft-core pair functions; needs no GPU.
add_executable(gapsys_pair_test gapsys_pair_test.cpp)

target_include_directories(gapsys_pair_test
  PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${CUDA_TOOLKIT_INCLUDE_DIR}"
    "${PROJECT_SOURCE_DIR}/cuda/include"
    "${PROJECT_SOURCE_DIR}/common/include"
    "${PROJECT_SOURCE_DIR}/cpu/include"
)

# Test the precision the core is built with.
if(QDYN_SPFP)
  target_compile_definitions(gapsys_pair_test PRIVATE QDYN_SPFP)
endif()

add_test(NAME gapsys_pair_test COMMAND gapsys_pair_test)
