#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "context.h"
#include "cpu_handler.h"
#include "cuda_handler.cuh"
#include "init.h"
#include "q_input.h"
#include "q_output.h"

static bool ends_with(const char *value, const char *suffix) {
    size_t value_len = strlen(value);
    size_t suffix_len = strlen(suffix);
    if (suffix_len > value_len) return false;
    return strcmp(value + value_len - suffix_len, suffix) == 0;
}

static void make_default_workdir(const char *input_file, char *workdir, size_t workdir_size) {
    char cwd[1024];
    char input_copy[1024];
    char *base;
    char stem[256];

    if (getcwd(cwd, sizeof(cwd)) == NULL) {
        printf(">>> FATAL: Could not determine current working directory. Exiting...\n");
        exit(EXIT_FAILURE);
    }

    strncpy(input_copy, input_file, sizeof(input_copy) - 1);
    input_copy[sizeof(input_copy) - 1] = '\0';
    base = basename(input_copy);
    strncpy(stem, base, sizeof(stem) - 1);
    stem[sizeof(stem) - 1] = '\0';
    char *dot = strrchr(stem, '.');
    if (dot != NULL) *dot = '\0';

    snprintf(workdir, workdir_size, "%s/.qgpu/%s", cwd, stem);
}

static void usage() {
    printf("qdyn [--gpu] input.inp\n");
    printf("qdyn [--gpu] --csv-dir input_folder\n");
}

static void calc_integration() {
    auto& ctx = Context::instance();
    ctx.init();
    if (q_input_mode == Q_INPUT_NATIVE) {
        q_output_init();
    } else {
        write_headers();
    }

    Handler& handler = ctx.run_gpu ? static_cast<Handler&>(CudaHandler::instance()) : static_cast<Handler&>(CpuHandler::instance());
    handler.initialize();
    handler.run(ctx.md.steps);
    q_output_finish();
    q_output_shutdown();
    handler.shutdown();
}

int main(int argc, char *argv[]) {
    const char *input_file = NULL;
    const char *csv_dir = NULL;
    bool requested_gpu = false;

    if (argc < 2) {
        usage();
        printf(">>> FATAL: Input .inp file or --csv-dir expected. Exiting...\n");
        return EXIT_FAILURE;
    }

    Context ctx;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--gpu") == 0) {
            requested_gpu = true;
        }
        else if (strcmp(argv[i], "--csv-dir") == 0) {
            if (i + 1 >= argc) {
                usage();
                printf(">>> FATAL: --csv-dir requires a folder argument. Exiting...\n");
                return EXIT_FAILURE;
            }
            csv_dir = argv[++i];
        }
        else if (input_file == NULL) {
            input_file = argv[i];
        }
        else {
            usage();
            printf(">>> FATAL: Unexpected argument '%s'. Exiting...\n", argv[i]);
            return EXIT_FAILURE;
        }
    }

    ctx.run_gpu = requested_gpu;

    if (csv_dir != NULL) {
        if (input_file != NULL) {
            usage();
            printf(">>> FATAL: --csv-dir cannot be combined with an .inp file. Exiting...\n");
            return EXIT_FAILURE;
        }
        q_input_mode = Q_INPUT_CSV;
        ctx.base_folder = csv_dir;
        calc_integration();
        return EXIT_SUCCESS;
    }

    if (input_file == NULL || !ends_with(input_file, ".inp")) {
        usage();
        printf(">>> FATAL: Q input mode requires an .inp file. Exiting...\n");
        return EXIT_FAILURE;
    }

    char workdir[1024];
    q_input_mode = Q_INPUT_NATIVE;
    make_default_workdir(input_file, workdir, sizeof(workdir));
    q_prepare_input(input_file, workdir);
    ctx.base_folder = workdir;
    calc_integration();

    return EXIT_SUCCESS;
}
