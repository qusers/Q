#pragma once

struct LincsSettings {
    int expansion_order = 8;
    int minimum_rotation_iterations = 2;
    int maximum_rotation_iterations = 8;
    double accuracy_tolerance = 1e-6;
};
