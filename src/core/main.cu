#include <stdio.h>
#include <stdlib.h>

#include <memory>

#include "command_parser.h"
#include "context.h"
#include "cpu_handler.h"
#include "cuda_handler.cuh"

void calc_integration() {
    auto& ctx = Context::instance();
    std::unique_ptr<Handler> handler =
        ctx.command_info.requested_gpu
            ? std::unique_ptr<Handler>(std::make_unique<CudaHandler>())
            : std::unique_ptr<Handler>(std::make_unique<CpuHandler>());
    handler->initialize();
    handler->run(ctx.md.steps);
    handler->shutdown();
}

int main(int argc, char* argv[]) {
    CommandParseResult parsed = parse_command(argc, argv);
    if (!parsed.ok) {
        print_command_usage();
        printf(">>> FATAL: %s Exiting...\n", parsed.error.c_str());
        return EXIT_FAILURE;
    }
    Context ctx;
    const CommandInfo& command = parsed.command;
    ctx.command_info = command;
    calc_integration();
    return EXIT_SUCCESS;
}
