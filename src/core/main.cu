#include <stdio.h>
#include <stdlib.h>

#include <string>

#include "command_parser.h"
#include "context.h"
#include "cpu_handler.h"
#include "cuda_handler.cuh"
#include "helpers.h"
#include "init.h"
#include "q_input.h"
#include "q_output.h"

static void calc_integration() {
    auto& ctx = Context::instance();
    ctx.init();
    if (q_input_mode == Q_INPUT_NATIVE) {
        q_output_init();
    } else {
        write_headers();
    }

    Handler& handler = ctx.command_info.requested_gpu ? static_cast<Handler&>(CudaHandler::instance()) : static_cast<Handler&>(CpuHandler::instance());
    handler.initialize();
    handler.run(ctx.md.steps);
    q_output_finish();
    q_output_shutdown();
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
