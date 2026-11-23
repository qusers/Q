#pragma once
#include <stdexcept>
#include <string>

struct QgpuError : std::runtime_error {
    using std::runtime_error::runtime_error;
};
struct InputError : QgpuError {
    using QgpuError::QgpuError;
};
struct ConstraintError : QgpuError {
    using QgpuError::QgpuError;
};
struct CudaError : QgpuError {
    using QgpuError::QgpuError;
};