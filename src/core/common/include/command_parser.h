#pragma once

#include <stdio.h>
#include <str_helpers.h>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

enum class CommandInputMode {
    NativeInp
};

struct CommandInfo {
    bool requested_gpu = false;
    CommandInputMode input_mode = CommandInputMode::NativeInp;
    std::vector<std::string> input_files;
    int n_replicates() const {
        return static_cast<int>(input_files.size());
    }
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
    std::string primary_input;
    std::vector<std::string> additional_inputs;
    int requested_replicates = 1;
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
                std::size_t parsed_characters = 0;
                requested_replicates = std::stoi(value, &parsed_characters);
                if (parsed_characters != value.size()) {
                    throw std::invalid_argument("trailing characters");
                }
            } catch (...) {
                result.error = "--replicates requires a positive integer argument.";
                return result;
            }

            if (requested_replicates <= 0) {
                result.error = "--replicates requires a positive integer argument.";
                return result;
            }
            if (requested_replicates > 65535) {
                result.error = "--replicates exceeds the CUDA grid y-dimension limit (65535).";
                return result;
            }
            replicate_count_was_set = true;
        } else if (arg == "--replica-input") {
            if (i + 1 >= argc) {
                result.error = "--replica-input requires an .inp file argument.";
                return result;
            }
            additional_inputs.push_back(argv[++i]);
        } else if (primary_input.empty()) {
            primary_input = arg;
        } else {
            result.error = "Unexpected argument '" + arg + "'.";
            return result;
        }
    }

    if (primary_input.empty() || !ends_with(primary_input, ".inp")) {
        result.error = "Q input mode requires an .inp file.";
        return result;
    }
    command.input_mode = CommandInputMode::NativeInp;

    if (replicate_count_was_set && !additional_inputs.empty()) {
        result.error = "--replicates cannot be combined with --replica-input.";
        return result;
    }

    for (const auto& input : additional_inputs) {
        if (!ends_with(input, ".inp")) {
            result.error = "Every --replica-input value must end with .inp.";
            return result;
        }
    }

    if (additional_inputs.empty()) {
        command.input_files.assign(
            static_cast<std::size_t>(requested_replicates),
            primary_input);
    } else {
        command.input_files.reserve(1 + additional_inputs.size());
        command.input_files.push_back(primary_input);
        command.input_files.insert(
            command.input_files.end(),
            additional_inputs.begin(),
            additional_inputs.end());
    }

    if (command.input_files.size() > 65535) {
        result.error = "The replica count exceeds the CUDA grid y-dimension limit (65535).";
        return result;
    }

    if (command.n_replicates() > 1 && !command.requested_gpu) {
        result.error = "Batched replicates currently require --gpu.";
        return result;
    }

    result.ok = true;
    result.command = command;
    return result;
}
