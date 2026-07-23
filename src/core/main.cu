#include <stdio.h>
#include <stdlib.h>

#include <memory>

#include "command_parser.h"
#include "context.h"
#include "cpu_handler.h"
#include "cuda_handler.cuh"

void calc_integration(const CommandParseResult& parsed) {
    std::unique_ptr<Handler> handler = parsed.command.requested_gpu
                                           ? std::unique_ptr<Handler>(std::make_unique<CudaHandler>())
                                           : std::unique_ptr<Handler>(std::make_unique<CpuHandler>());
    handler->initialize(parsed.command);
    handler->run();
}

int main(int argc, char* argv[]) {
    CommandParseResult parsed = parse_command(argc, argv);
    if (!parsed.ok) {
        print_command_usage();
        printf(">>> FATAL: %s Exiting...\n", parsed.error.c_str());
        return EXIT_FAILURE;
    }
    if (parsed.command.n_replicates > 1) {
        printf(">>> FATAL: Multi-replica runtime support is not implemented yet. Exiting...\n");
        return EXIT_FAILURE;
    }
    calc_integration(parsed);
    return EXIT_SUCCESS;
}
