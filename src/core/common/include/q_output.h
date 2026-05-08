#ifndef __Q_OUTPUT_H__
#define __Q_OUTPUT_H__

#include "md_types.h"

void q_output_configure(const char *final_file,
                        const char *trajectory_file,
                        const char *energy_file,
                        const char *trajectory_atoms,
                        const char *topology_file,
                        int steps,
                        int trajectory_interval,
                        int energy_interval);

void q_output_init();
void q_output_step(int iteration);
void q_output_finish();
void q_output_shutdown();
void q_output_write_restart(const coord_t *final_coords, const vel_t *final_velocities, int natoms);

#endif /* __Q_OUTPUT_H__ */
