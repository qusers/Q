#ifndef __PARSE_H__
#define __PARSE_H__

#define N_COLUMNS       10
#define COLUMN_WIDTH    25

struct csvfile_t {
    char ***buffer;
    int n_lines;
    int ext;
};

csvfile_t read_csv(const char *filename, int ext, char *base_folder);
void clean_csv(csvfile_t file);

/* =============================================
 * == FROM MD FILE
 * =============================================
 */

void init_md(const char *filename);

/* =============================================
 * == FROM TOPOLOGY FILE
 * =============================================
 */

void init_topo(const char *filename);

void init_coords(const char *filename);
void init_bonds(const char *filename);
void init_cbonds(const char *filename);
void init_angles(const char *filename);
void init_cangles(const char *filename);
void init_torsions(const char *filename);
void init_ctorsions(const char *filename);
void init_impropers(const char *filename);
void init_cimpropers(const char *filename);
void init_charges(const char *filename);
void init_ccharges(const char *filename);
void init_LJ_matrix();
void init_ngbrs14(const char *filename);
void init_ngbrs23(const char *filename);
void init_ngbrs14_long(const char* filename);
void init_ngbrs23_long(const char* filename);
void init_catypes(const char *filename);
void init_atypes(const char *filename);
void init_excluded(const char *filename);
void init_molecules(const char *filename);
void init_charge_groups(const char *filename);

/* =============================================
 * == FROM FEP FILE
 * =============================================
 */

void init_qangcouples(const char *filename);
void init_qatoms(const char *filename);
void init_qcangles(const char *filename);
void init_qcatypes(const char *filename);
void init_qcbonds(const char *filename);
void init_qcimpropers(const char *filename);
void init_qctorsions(const char *filename);
void init_qoffdiags(const char *filename);
void init_qimprcouples(const char *filename);
void init_qsoftpairs(const char *filename);
void init_qtorcouples(const char *filename);

void init_qangles(const char *filename);
void init_qatypes(const char *filename);
void init_qbonds(const char *filename);
void init_qcharges(const char *filename);
void init_qelscales(const char *filename);
void init_qexclpairs(const char *filename);
void init_qimpropers(const char *filename);
void init_qshakes(const char *filename);
void init_qsoftcores(const char *filename);
void init_qtorsions(const char *filename);

/* =============================================
 * == FROM INPUT FILE
 * =============================================
 */

void init_icoords(const char *filename);
void init_ivelocities(const char *filename);    

#endif /* __PARSE_H__ */