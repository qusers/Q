#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <filesystem>
#include <string>
#include <system_error>

inline std::string make_default_workdir(const std::string& input_file) {
    std::error_code error;
    std::filesystem::path cwd = std::filesystem::current_path(error);
    if (error) {
        printf(">>> FATAL: Could not determine current working directory. Exiting...\n");
        exit(EXIT_FAILURE);
    }

    std::filesystem::path input_path(input_file);
    return (cwd / ".qgpu" / input_path.stem()).string();
}

static void fatal(const std::string& message) {
    printf(">>> FATAL: %s\n", message.c_str());
    exit(EXIT_FAILURE);
}

static bool file_exists(const std::string& path) {
    return access(path.c_str(), F_OK) == 0;
}

static std::string dirname_of(const std::string& path) {
    size_t slash = path.find_last_of('/');
    if (slash == std::string::npos) return ".";
    if (slash == 0) return "/";
    return path.substr(0, slash);
}

static bool is_absolute(const std::string& path) {
    return !path.empty() && path[0] == '/';
}

static std::string join_path(const std::string& left, const std::string& right) {
    if (left.empty() || left == ".") return right;
    if (!right.empty() && right[0] == '/') return right;
    if (left[left.size() - 1] == '/') return left + right;
    return left + "/" + right;
}

static std::string resolve_existing_or_cwd(const std::string& value, const std::string& input_dir) {
    if (value.empty()) return "";
    if (is_absolute(value)) return value;
    std::string input_relative = join_path(input_dir, value);
    if (file_exists(input_relative)) return input_relative;
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) == NULL) fatal("Could not determine current working directory.");
    return join_path(cwd, value);
}

static std::string resolve_output_path(const std::string& value, const std::string& input_dir) {
    if (value.empty()) return "";
    if (is_absolute(value)) return value;
    return join_path(input_dir, value);
}