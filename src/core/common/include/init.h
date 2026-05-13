#pragma once

#define __PROFILING__
// #define DEBUG
#define VERBOSE

void exclude_qatom_definitions();
void init_pshells_from_charge_groups();
void exclude_shaken_definitions();
void init_unified_atom_parameters();
void init_patoms();
void initial_shaking();
void stop_cm_translation();
void finalize_ngbrs14();
void init_atoms_list();
void initialize_charge_tables();
void initialize_catype_tables();
void init_inv_mass();
void init_pshells();
void init_shake();
void init_velocities();
void init_water_sphere();
void init_wshells();

