#ifndef __PARSE_H__
#define __PARSE_H__

#define N_COLUMNS       10
#define COLUMN_WIDTH    25

struct csvfile_t {
    char ***buffer;
    int n_lines;
    int ext;
};

csvfile_t read_csv(char *filename, int ext, char *base_folder);
void clean_csv(csvfile_t file);

/* =============================================
 * == FROM MD FILE
 * =============================================
 */

void init_md(char *filename);

/* =============================================
 * == FROM TOPOLOGY FILE
 * =============================================
 */

void init_topo(char *filename);

void init_coords(char* filename);
void init_bonds(char* filename);
void init_cbonds(char* filename);
void init_angles(char* filename);
void init_cangles(char* filename);
void init_torsions(char* filename);
void init_ctorsions(char* filename);
void init_impropers(char* filename);
void init_cimpropers(char* filename);
void init_charges(char* filename);
void init_ccharges(char* filename);
void init_ngbrs14(char* filename);
void init_ngbrs23(char* filename);
void init_catypes(char* filename);
void init_atypes(char* filename);
void init_excluded(char* filename);

#endif /* __PARSE_H__ */