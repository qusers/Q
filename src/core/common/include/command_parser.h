#pragma once

#include <stdio.h>
#include <str_helpers.h>

#include <string>

struct CommandInfo {
    bool requested_gpu = false;
    std::string input_file;
};

struct CommandParseResult {
    bool ok = false;
    CommandInfo command;
    std::string error;
};

inline void print_command_usage() {
    printf("qdyn [--gpu] input.inp\n");
}

inline CommandParseResult parse_command(int argc, char* argv[]) {
    CommandParseResult result;

    if (argc < 2) {
        result.error = "Input .inp file expected.";
        return result;
    }

    CommandInfo command;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "--gpu") {
            command.requested_gpu = true;
        } else if (command.input_file.empty()) {
            command.input_file = arg;
        } else {
            result.error = "Unexpected argument '" + arg + "'.";
            return result;
        }
    }

    if (command.input_file.empty() || !ends_with(command.input_file, ".inp")) {
        result.error = "Q input requires an .inp file.";
        return result;
    }

    result.ok = true;
    result.command = command;
    return result;
}
