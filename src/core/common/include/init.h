#pragma once
#include "shake.h"
class Context;
void exclude_qatom_definitions();
void init_pshells_from_charge_groups();
void exclude_shaken_definitions(Context& ctx);
void init_unified_atom_parameters();
void init_patoms();
void finalize_ngbrs14();
void init_atoms_list();
void initialize_charge_tables();
void initialize_catype_tables();
void init_inv_mass();
void init_pshells();
void init_shake(Context& ctx);
void init_velocities();
void init_water_sphere();
void init_wshells();

void stop_cm_translation(Context& ctx);
void init_for_temperature(Context& ctx, Shake &shake);