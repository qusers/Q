#include <stdio.h>
#include <stdlib.h>

#include <memory>

#include "command_parser.h"
#include "context.h"
#include "run_simulation.h"

int main(int argc, char* argv[]) {
    CommandParseResult parsed = parse_command(argc, argv);
    if (!parsed.ok) {
        print_command_usage();
        printf(">>> FATAL: %s Exiting...\n", parsed.error.c_str());
        return EXIT_FAILURE;
    }

    try {
        run_simulation(parsed.command);
    } catch (const std::exception& e) {
        printf(">>> FATAL: %s Exiting...\n", e.what());
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
