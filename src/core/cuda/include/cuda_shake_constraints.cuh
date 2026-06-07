#pragma once

class Context;
void init_shake_constraints_kernel_data(Context& ctx);

void calc_shake_constraints_host(Context& ctx);
void cleanup_shake_constraints();
