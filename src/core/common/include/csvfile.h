struct csvfile_t {
    char*** buffer = nullptr;
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

csvfile_t read_csv(const char* filename, int ext, const char* base_folder);