#include "system.h"
#include "nonbonded.h"

/* =============================================
 * == NON-BONDED INTERACTIONS
 * =============================================
 */

 void calc_nonbonded_pp_forces() {
    bool bond14, bond23;
    double scaling;
    coord_t da;
    double r2a, ra, r6a;
    double Vela, V_a, V_b;
    double dva;
    double qi, qj;
    double ai_aii, aj_aii, ai_bii, aj_bii;
    int i, j;
    catype_t ai_type, aj_type;

    if (n_patoms > 0) {
        energies.Ucoul = 0;
        energies.Uvdw = 0;
    }

    for (int pi = 0; pi < n_patoms; pi++) {
        for (int pj = i+1; pj < n_patoms; pj++) {
            i = p_atoms[pi].a - 1;
            j = p_atoms[pj].a - 1;
            bond23 = false;
            bond14 = false;
            for (int k = 0; k < n_ngbrs23; k++) {
                if (ngbrs23[k].ai == charges[i].a && ngbrs23[k].aj == charges[j].a) {
                    // printf("BOND 23: %d, %d", ngbrs23[k].ai, ngbrs23[k].aj);
                    bond23 = true;
                }
            }
            for (int k = 0; k < n_ngbrs14; k++) {
                if (ngbrs14[k].ai == charges[i].a && ngbrs14[k].aj == charges[j].a) {
                    // printf("BOND 14: %d, %d", ngbrs23[k].ai, ngbrs23[k].aj);
                    bond14 = true;
                }
            }

            if (bond23) continue;

            scaling = bond14 ? .5 : 1;

            qi = ccharges[charges[i].code - 1].charge;
            qj = ccharges[charges[j].code - 1].charge;

            ai_type = catypes[atypes[i].code - 1];
            aj_type = catypes[atypes[j].code - 1];

            da.x = coords[j].x - coords[i].x;
            da.y = coords[j].y - coords[i].y;
            da.z = coords[j].z - coords[i].z;
            r2a = 1 / (pow(da.x, 2) + pow(da.y, 2) + pow(da.z, 2));
            ra = sqrt(r2a);
            r6a = r2a * r2a * r2a;

            Vela = scaling * Coul * qi * qj * ra;

            ai_aii = bond14 ? ai_type.aii_1_4 : ai_type.aii_normal;
            aj_aii = bond14 ? aj_type.aii_1_4 : aj_type.aii_normal;
            ai_bii = bond14 ? ai_type.bii_1_4 : ai_type.bii_normal;
            aj_bii = bond14 ? aj_type.bii_1_4 : aj_type.bii_normal;

            V_a = r6a * r6a * ai_aii * aj_aii;
            V_b = r6a * ai_bii * aj_bii;
            dva = r2a * ( -Vela -12 * V_a + 6 * V_b);

            dvelocities[i].x -= dva * da.x;
            dvelocities[i].y -= dva * da.y;
            dvelocities[i].z -= dva * da.z;

            dvelocities[j].x += dva * da.x;
            dvelocities[j].y += dva * da.y;
            dvelocities[j].z += dva * da.z;

            energies.Ucoul += Vela;
            energies.Uvdw += (V_a - V_b);
        }
    }
}
