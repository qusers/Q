#ifndef __PARSE_H__
#define __PARSE_H__

#include <vector>

struct csvfile_t {
    char ***buffer = nullptr;
    int n_lines = 0;
    int ext = 0;

    csvfile_t() = default;
    csvfile_t(const csvfile_t&) = delete;
    csvfile_t& operator=(const csvfile_t&) = delete;
    csvfile_t(csvfile_t&& other) noexcept;
    csvfile_t& operator=(csvfile_t&& other) noexcept;
    ~csvfile_t();

    void reset() noexcept;
};

struct charge_group_t {
    int iswitch = 0;
    std::vector<int> atoms;
};

struct charge_group_config_t {
    int n_cgrps_solute = 0;
    int n_cgrps_solvent = 0;
    int iuse_switch_atom = 0;
    std::vector<charge_group_t> charge_groups;
};

csvfile_t read_csv(const char *filename, int ext, const char *base_folder);

/* =============================================
 * == FROM MD FILE
 * =============================================
 */

void parse_md(const char *filename);

/* =============================================
 * == FROM TOPOLOGY FILE
 * =============================================
 */

void parse_topo(const char *filename);

void parse_coords(const char *filename);
void parse_bonds(const char *filename);
void parse_cbonds(const char *filename);
void parse_angles(const char *filename);
void parse_cangles(const char *filename);
void parse_torsions(const char *filename);
void parse_ctorsions(const char *filename);
void parse_impropers(const char *filename);
void parse_cimpropers(const char *filename);
void parse_charges(const char *filename);
void parse_ccharges(const char *filename);
void parse_LJ_matrix();
void parse_ngbrs14(const char *filename);
void parse_ngbrs23(const char *filename);
void parse_ngbrs14_long(const char* filename);
void parse_ngbrs23_long(const char* filename);
void parse_catypes(const char *filename);
void parse_atypes(const char *filename);
void parse_excluded(const char *filename);
void parse_molecules(const char *filename);
charge_group_config_t read_charge_groups(const char *filename);

/* =============================================
 * == FROM FEP FILE
 * =============================================
 */

void parse_qangcouples(const char *filename);
void parse_qatoms(const char *filename);
void parse_qcangles(const char *filename);
void parse_qcatypes(const char *filename);
void parse_qcbonds(const char *filename);
void parse_qcimpropers(const char *filename);
void parse_qctorsions(const char *filename);
void parse_qoffdiags(const char *filename);
void parse_qimprcouples(const char *filename);
void parse_qsoftpairs(const char *filename);
void parse_qtorcouples(const char *filename);

void parse_qangles(const char *filename);
void parse_qatypes(const char *filename);
void parse_qbonds(const char *filename);
void parse_qcharges(const char *filename);
void parse_qelscales(const char *filename);
void parse_qexclpairs(const char *filename);
void parse_qimpropers(const char *filename);
void parse_qshakes(const char *filename);
void parse_qsoftcores(const char *filename);
void parse_qtorsions(const char *filename);

/* =============================================
 * == FROM INPUT FILE
 * =============================================
 */

void parse_icoords(const char *filename);
void parse_ivelocities(const char *filename);    

#endif /* __PARSE_H__ */
