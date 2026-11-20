#pragma once
#include <vector>

#include "parse.h"

struct StageSettings {
    md_t md;                      // only "Stage" fields are applied
    std::vector<double> lambdas;  // size must equal number of FEP states
    std::vector<restrpos_t> restrspos;
    std::vector<restrseq_t> restrseqs;
    std::vector<restrdis_t> restrdists;
    std::vector<restrang_t> restrangs;
    std::vector<restrwall_t> restrwalls;
    NativeOutputConfig outputs;  // empty path = don't write that file
};

struct SystemState {
    std::vector<coord_t> coords;
    std::vector<vel_t> velocities;
    std::vector<double> theta_corr;  // water polarisation shells
};

StageSettings stage_from(const ParseResult& p);  // copy md, lambdas, restr*, native_output
SystemState state_from(const ParseResult& p);    // copy coords, velocities, restart_theta_corr