#pragma once

#include <stdio.h>
#include <str_helpers.h>

#include <stdexcept>
#include <string>
#include <vector>

enum class CommandInputMode {
    NativeInp
};

struct CommandInfo {
    bool requested_gpu = false;
    std::string input_file;
    CommandInputMode input_mode = CommandInputMode::NativeInp;
    int n_replicates = 1;
    std::vector<std::string> replica_input_files;
};

struct CommandParseResult {
    bool ok = false;
    CommandInfo command;
    std::string error;
};

inline void print_command_usage() {
    printf("qdyn [--gpu] input.inp\n");
    printf("qdyn --gpu --replicates N input.inp\n");
    printf("qdyn --gpu input0.inp [--replica-input input1.inp ...]\n");
}

inline CommandParseResult parse_command(int argc, char* argv[]) {
    CommandParseResult result;

    if (argc < 2) {
        result.error = "Input .inp file expected.";
        return result;
    }

    CommandInfo command;
    bool replicate_count_was_set = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "--gpu") {
            command.requested_gpu = true;
        } else if (arg == "--replicates") {
            if (replicate_count_was_set) {
                result.error = "--replicates may only be specified once.";
                return result;
            }
            if (i + 1 >= argc) {
                result.error = "--replicates requires a positive integer argument.";
                return result;
            }

            try {
                std::string value = argv[++i];
                size_t parsed_characters = 0;
                command.n_replicates = std::stoi(value, &parsed_characters);
                if (parsed_characters != value.size()) {
                    throw std::invalid_argument("trailing characters");
                }
            } catch (...) {
                result.error = "--replicates requires a positive integer argument.";
                return result;
            }

            if (command.n_replicates <= 0) {
                result.error = "--replicates requires a positive integer argument.";
                return result;
            }
            if (command.n_replicates > 65535) {
                result.error = "--replicates exceeds the CUDA grid y-dimension limit (65535).";
                return result;
            }
            replicate_count_was_set = true;
        } else if (arg == "--replica-input") {
            if (i + 1 >= argc) {
                result.error = "--replica-input requires an .inp file argument.";
                return result;
            }
            command.replica_input_files.push_back(argv[++i]);
        } else if (command.input_file.empty()) {
            command.input_file = arg;
        } else {
            result.error = "Unexpected argument '" + arg + "'.";
            return result;
        }
    }

    if (command.input_file.empty() || !ends_with(command.input_file, ".inp")) {
        result.error = "Q input mode requires an .inp file.";
        return result;
    }
    command.input_mode = CommandInputMode::NativeInp;

    if (!command.replica_input_files.empty()) {
        if (replicate_count_was_set) {
            result.error = "--replicates cannot be combined with --replica-input.";
            return result;
        }
        for (const auto& input : command.replica_input_files) {
            if (!ends_with(input, ".inp")) {
                result.error = "Every --replica-input value must end with .inp.";
                return result;
            }
        }
        command.n_replicates =
            1 + static_cast<int>(command.replica_input_files.size());
        if (command.n_replicates > 65535) {
            result.error = "The replica count exceeds the CUDA grid y-dimension limit (65535).";
            return result;
        }
    }

    if (command.n_replicates > 1 && !command.requested_gpu) {
        result.error = "Batched replicates currently require --gpu.";
        return result;
    }

    result.ok = true;
    result.command = command;
    return result;
}
