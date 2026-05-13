#pragma once

static void fatal(const std::string& message) {
    printf(">>> FATAL: %s\n", message.c_str());
    exit(EXIT_FAILURE);
}
