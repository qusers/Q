#pragma once
#include <string>

#include "qgpu_error.h"

[[noreturn]] inline void fatal(const std::string& message) {
    throw QgpuError(message);
}
