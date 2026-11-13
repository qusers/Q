#pragma once

#include <memory>

#include "command_parser.h"
#include "cpu_handler.h"
#include "cuda_handler.cuh"

// Shared by the CLI and Python extension.
// GPU selection, when needed, happens before this call.
inline void run_simulation(const CommandInfo& command) {
    std::unique_ptr<Handler> handler;

    if (command.requested_gpu) {
        handler = std::make_unique<CudaHandler>();
    } else {
        handler = std::make_unique<CpuHandler>();
    }

    handler->initialize(command);
    handler->run();
}