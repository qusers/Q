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

#endif /* __PARSE_H__ */