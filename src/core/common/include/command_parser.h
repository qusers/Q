#pragma once

#include <stdio.h>
#include <str_helpers.h>

#include <string>

enum class CommandInputMode {
    Csv,
    NativeInp
};

struct CommandInfo {
    bool requested_gpu = false;
    std::string input_file;
    std::string csv_dir;
    CommandInputMode input_mode = CommandInputMode::NativeInp;
};

struct CommandParseResult {
    bool ok = false;
    CommandInfo command;
    std::string error;
};

inline void print_command_usage() {
    printf("qdyn [--gpu] input.inp\n");
    printf("qdyn [--gpu] --csv-dir input_folder\n");
}

inline CommandParseResult parse_command(int argc, char* argv[]) {
    CommandParseResult result;

    if (argc < 2) {
        result.error = "Input .inp file or --csv-dir expected.";
        return result;
    }

    CommandInfo command;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "--gpu") {
            command.requested_gpu = true;
        } else if (arg == "--csv-dir") {
            if (i + 1 >= argc) {
                result.error = "--csv-dir requires a folder argument.";
                return result;
            }
            command.csv_dir = argv[++i];
        } else if (command.input_file.empty()) {
            command.input_file = arg;
        } else {
            result.error = "Unexpected argument '" + arg + "'.";
            return result;
        }
    }

    if (!command.csv_dir.empty()) {
        if (!command.input_file.empty()) {
            result.error = "--csv-dir cannot be combined with an .inp file.";
            return result;
        }
        command.input_mode = CommandInputMode::Csv;
    } else {
        if (command.input_file.empty() || !ends_with(command.input_file, ".inp")) {
            result.error = "Q input mode requires an .inp file.";
            return result;
        }
        command.input_mode = CommandInputMode::NativeInp;
    }

    result.ok = true;
    result.command = command;
    return result;
}
