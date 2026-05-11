#pragma once
#include <algorithm>
#include <cctype>
#include <sstream>
#include <string>
#include <vector>

#include "helpers.h"

inline bool ends_with(const std::string& value, const std::string& suffix) {
    if (suffix.size() > value.size()) return false;
    return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static std::string trim(const std::string& value) {
    size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

static std::string strip_comment(const std::string& line) {
    size_t bang = line.find('!');
    if (bang == std::string::npos) return trim(line);
    return trim(line.substr(0, bang));
}

static std::string lower_normalized(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        if (c == '_') return '-';
        return (char)std::tolower(c);
    });
    return value;
}

static std::vector<std::string> split_ws(const std::string& line) {
    std::vector<std::string> fields;
    std::istringstream iss(line);
    std::string field;
    while (iss >> field) fields.push_back(field);
    return fields;
}

static std::vector<std::vector<std::string> > split_groups(const std::vector<std::string>& flat, size_t width) {
    if (width == 0) fatal("Invalid grouped input width.");
    if (flat.size() % width != 0) {
        fatal("Invalid grouped input length: " + std::to_string(flat.size()) + " values cannot be split into groups of " + std::to_string(width) + ".");
    }

    std::vector<std::vector<std::string> > groups;
    for (size_t i = 0; i + width <= flat.size(); i += width) {
        groups.push_back(std::vector<std::string>(flat.begin() + i, flat.begin() + i + width));
    }
    return groups;
}

std::string strip_cr(std::string value) {
    if (!value.empty() && value.back() == '\r') {
        value.pop_back();
    }
    return value;
}

std::vector<std::string> split_semicolon(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    for (char ch : line) {
        if (ch == ';') {
            fields.push_back(strip_cr(field));
            field.clear();
        } else {
            field.push_back(ch);
        }
    }
    fields.push_back(strip_cr(field));
    return fields;
}

int to_int(const std::string& value) {
    return std::atoi(value.c_str());
}

double to_double(const std::string& value) {
    return std::strtod(value.c_str(), nullptr);
}