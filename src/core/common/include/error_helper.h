#pragma once

#include <stdexcept>

std::runtime_error parse_error(const std::string& message) {
    return std::runtime_error(message);
}