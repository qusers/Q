#include <stdio.h>
#include <stdlib.h>

#include "command_parser.h"
#include "context.h"
#include "cpu_handler.h"
#include "cuda_handler.cuh"
#include "debug.h"

static void calc_integration() {
    auto& ctx = Context::instance();
    ctx.init();
    // output_ctx_in_file(ctx);
    Handler& handler = ctx.command_info.requested_gpu ? static_cast<Handler&>(CudaHandler::instance()) : static_cast<Handler&>(CpuHandler::instance());
    handler.initialize();
    handler.run(ctx.md.steps);
    handler.shutdown();
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
