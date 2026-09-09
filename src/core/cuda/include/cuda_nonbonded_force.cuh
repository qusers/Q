#pragma once
#include "nonbonded_force.h"
#include "precision.h"

enum GroupPairMode : uint8_t {
    GROUP_PAIR_IGNORE = 0,
    GROUP_PAIR_EXACT = 1,
    GROUP_PAIR_LRF = 2
};

struct ExactEntry {
    int x_start;
    int x_len;  // <= 32
    int y_start;
    int y_len;  // <= 32

    uint8_t diagonal;
    uint8_t y_indirect;
    int y_type_slot;
};

struct LrfPairEntry {
    int range1;
    int range2;

    int group1;
    int group2;
};

class CudaNonbondedForce final : public NonbondedForce {
   public:
    void calc(Context& ctx) override;

   protected:
    void init_backend(Context& ctx) override;

   private:
    void calc_all_direct_pairs(Context& ctx);
    void init_calculation_groups(Context& ctx);
    void init_calculation_groups_by_switch(Context& ctx);
    void init_calculation_groups_by_all_atoms(Context& ctx);
    void init_lrf_coefficients(Context& ctx);
    void calc_exact_tiles(Context& ctx);
    void calc_lrf(Context& ctx);
    void build_lrf_atom_csr(Context& ctx);
    void build_exact_atom_tiles(Context& ctx);

    std::unique_ptr<HostDeviceBuffer<real_t>> coord_x_, coord_y_, coord_z_;

    std::unique_ptr<HostDeviceBuffer<uint8_t>> group_pair_modes_;

    std::unique_ptr<HostDeviceBuffer<ExactEntry>> exact_tiles_;
    std::unique_ptr<HostDeviceBuffer<LrfPairEntry>> lrf_group_pairs_;

    std::unique_ptr<HostDeviceBuffer<int>> exact_tile_count_;
    std::unique_ptr<HostDeviceBuffer<int>> lrf_pair_count_;
    std::unique_ptr<HostDeviceBuffer<int>> list_overflow_;

    std::unique_ptr<HostDeviceBuffer<LrfCoefficients>> lrf_coefficients_;

    size_t exact_tile_capacity_ = 0;
    size_t lrf_pair_capacity_ = 0;

    int n_exact_tiles_ = 0;
    int n_lrf_pairs_ = 0;

    std::unique_ptr<HostDeviceBuffer<int>> lrf_atom_degrees_;

    std::unique_ptr<HostDeviceBuffer<int>> lrf_atom_offsets_;

    std::unique_ptr<HostDeviceBuffer<int>> lrf_source_atom_slots_;

    std::unique_ptr<HostDeviceBuffer<unsigned char>> lrf_scan_temp_;

    size_t lrf_source_atom_capacity_ = 0;
    size_t lrf_scan_temp_bytes_ = 0;
    int n_lrf_source_atom_entries_ = 0;

    /*
     * Geometry associated with every CSR source-atom entry.
     */
    std::unique_ptr<HostDeviceBuffer<double>> lrf_dx_;
    std::unique_ptr<HostDeviceBuffer<double>> lrf_dy_;
    std::unique_ptr<HostDeviceBuffer<double>> lrf_dz_;

    std::unique_ptr<HostDeviceBuffer<double>> lrf_q_r1_;
    std::unique_ptr<HostDeviceBuffer<double>> lrf_q_r3_;
    std::unique_ptr<HostDeviceBuffer<double>> lrf_q_r5_;
    std::unique_ptr<HostDeviceBuffer<double>> lrf_q_r7_;

    std::unique_ptr<HostDeviceBuffer<int>> exact_atom_degrees_;
    std::unique_ptr<HostDeviceBuffer<int>> exact_atom_offsets_;
    std::unique_ptr<HostDeviceBuffer<int>> exact_entry_degrees_;
    std::unique_ptr<HostDeviceBuffer<int>> exact_entry_offsets_;

    std::unique_ptr<HostDeviceBuffer<int>> exact_source_atom_slots_;
    size_t exact_source_atom_capacity_ = 0;
    int n_exact_source_atoms_ = 0;
};
