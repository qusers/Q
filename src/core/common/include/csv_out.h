#pragma once

#include <string>

#include "base_output.h"

class CsvOutput : public BaseOutput {
   public:
    explicit CsvOutput(const std::string& csv_dir);
    ~CsvOutput() override = default;

    void init(Context& ctx) override;

   protected:
    void output_trajectory(Context& ctx, int iteration) override;
    void output_energy(Context& ctx, int iteration) override;

   private:
    std::string csv_dir_;

    std::string output_path(const char* filename) const;
    void write_atom_count_header(Context& ctx, const char* filename) const;
    void write_energy_header(Context& ctx) const;
};
