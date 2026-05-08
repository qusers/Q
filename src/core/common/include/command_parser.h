#include <stdio.h>

#include <string>
#pragma once

struct CommandInfo {
    bool requested_gpu;
    std::string* input_files;
    std::string* csv_dir;
};

static void helper() {
    printf("qdyn [--gpu] input.inp\n");
    printf("qdyn [--gpu] --csv-dir input_folder\n");
}

void parse_command(int argc, char* argv[]) {
    
}