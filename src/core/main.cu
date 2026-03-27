#include "cpu/include/context.h"
#include <math.h>
#include <stdio.h>

#include "system.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf(">>> FATAL: Input file folder expected. Exiting...\n");
        exit(EXIT_FAILURE);
    }
    auto &ctx = Context::instance();
    bool run_gpu;

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
