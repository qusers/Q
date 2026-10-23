#pragma once

#include <cmath>

#include "context.h"
#include "geometry.h"

HD inline offdiag_record_t calc_offdiag(
    const offdiag_term_t& term, const coord_t* coords) {
    const coord_t d = coords[term.al] - coords[term.ak];
    const double r = norm(d);
    const double Hij = term.A * exp(-term.mu * r);

    return {term.state_i, term.state_j, Hij, r};
}

void calc_offdiags_cpu(Context& ctx);
void calc_offdiags_cuda(Context& ctx);