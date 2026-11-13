#pragma once

#include <stdexcept>

#include "qgpu_error.h"

inline InputError parse_error(const std::string& message) {
    return InputError(message);
}