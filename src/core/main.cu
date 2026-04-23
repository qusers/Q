#include <math.h>
#include <stdio.h>

#include "context.h"
#include "cpu_handler.h"
#include "cuda_handler.cuh"
#include "init.h"

void calc_integration() {
    auto& ctx = Context::instance();
    ctx.init();
    init_variables();
    Handler& handler = ctx.run_gpu ? static_cast<Handler&>(CudaHandler::instance()) : static_cast<Handler&>(CpuHandler::instance());
    handler.initialize();
    handler.run(ctx.md.steps);
    handler.shutdown();
    clean_variables();
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf(">>> FATAL: Input file folder expected. Exiting...\n");
        exit(EXIT_FAILURE);
    }
    auto& ctx = Context::instance();

    if (argc > 2) {
        ctx.run_gpu = strcmp(argv[1], "--gpu") == 0;
        ctx.base_folder = argv[2];
    } else {
        ctx.run_gpu = false;
        ctx.base_folder = argv[1];
    }

    calc_integration();

    exit(EXIT_SUCCESS);
}
