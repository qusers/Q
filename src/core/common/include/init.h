#pragma once

#define __PROFILING__
//#define DEBUG
#define VERBOSE


void init_variables();
void clean_variables();

/* =============================================
 * == DEVICE SETTINGS
 * =============================================
 */

/* =============================================
 * == GENERAL
 * =============================================
 */

/* =============================================
 * == FROM MD FILE
 * =============================================
 */

void init_pshells();
void init_pshells_with_switch_atoms();
void init_pshells_with_centroids();
void init_restrseqs(char* filename);


void init_water_sphere();
void init_wshells();

/* =============================================
 * == SHAKE
 * =============================================
 */

void init_shake();


void init_velocities();
void init_dvelocities();
void init_energies();

/* =============================================
 * == ENERGY & TEMPERATURE
 * =============================================
 */


/* =============================================
 * == INTEGRATION METHODS
 * =============================================
 */

void init_variables();
void clean_variables();
// void write_header(const char *filename);
// void write_energy_header();
// void write_coords(int iteration);
// void write_velocities(int iteration);
// void write_energies(int iteration);

void calc_integration();
void calc_integration_step(int iteration);

