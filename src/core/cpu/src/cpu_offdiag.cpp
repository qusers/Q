#include "offdiag.h"

void calc_offdiags_cpu(Context& ctx) {
    const int n = ctx.n_offdiags();
    if (n == 0) return;

    const auto* terms = ctx.offdiag_terms->cpu_data_p;
    const auto* coords = ctx.coords->cpu_data_p;
    auto* records = ctx.offdiag_records->cpu_data_p;

    for (int i = 0; i < n; ++i) {
        records[i] = calc_offdiag(terms[i], coords);
    }
}