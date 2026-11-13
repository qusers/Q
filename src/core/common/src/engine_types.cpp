#include "engine_types.h"

StageSettings stage_from(const ParseResult& p) {
    StageSettings s;
    s.md = p.md;
    s.lambdas = p.lambdas;
    s.restrspos = p.restrspos;
    s.restrseqs = p.restrseqs;
    s.restrdists = p.restrdists;
    s.restrangs = p.restrangs;
    s.restrwalls = p.restrwalls;
    s.outputs = p.native_output;
    return s;
}
SystemState state_from(const ParseResult& p) {
    return {p.coords, p.velocities, p.restart_theta_corr};
}