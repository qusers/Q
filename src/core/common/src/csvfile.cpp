#include "csvfile.h"

#include <stdio.h>
#include <unistd.h>

#include <cmath>
#include <cstring>

#include "constants.h"
csvfile_t::csvfile_t(csvfile_t&& other) noexcept
    : buffer(other.buffer), n_lines(other.n_lines), ext(other.ext) {
    other.buffer = nullptr;
    other.n_lines = 0;
    other.ext = 0;
}

csvfile_t& csvfile_t::operator=(csvfile_t&& other) noexcept {
    if (this != &other) {
        reset();
        buffer = other.buffer;
        n_lines = other.n_lines;
        ext = other.ext;

        other.buffer = nullptr;
        other.n_lines = 0;
        other.ext = 0;
    }

    return *this;
}

csvfile_t::~csvfile_t() {
    reset();
}

void csvfile_t::reset() noexcept {
    if (buffer == nullptr || n_lines <= 0) {
        buffer = nullptr;
        n_lines = 0;
        ext = 0;
        return;
    }

    for (int i = 0; i <= n_lines + ext; i++) {
        for (int j = 0; j < N_COLUMNS; j++) {
            free(buffer[i][j]);
        }
        free(buffer[i]);
    }
    free(buffer);

    buffer = nullptr;
    n_lines = 0;
    ext = 0;
}

csvfile_t read_csv(const char* filename, int ext, const char* base_folder) {
    csvfile_t retval;

    retval.ext = ext;

    char path[1024];
    sprintf(path, "%s/%s", base_folder, filename);
    if (access(path, F_OK) == -1) {
        printf(">>> FATAL: The following file could not be found. Exiting...\n");
        puts(path);
        exit(EXIT_FAILURE);
    }

    // File handle
    FILE* fp;

    fp = fopen(path, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    // Get number of lines
    char nlines[COLUMN_WIDTH];

    if (fgets(nlines, COLUMN_WIDTH, fp)) {
        retval.n_lines = atoi(nlines);
    } else {
        retval.n_lines = 0;
        return retval;
    }

    if (retval.n_lines == 0) {
        return retval;
    }

    char line[N_COLUMNS * COLUMN_WIDTH];
    retval.buffer = (char***)malloc(retval.n_lines * N_COLUMNS * sizeof(char**));

    for (int i = 0; i <= retval.n_lines + ext; i++) {
        retval.buffer[i] = (char**)malloc(N_COLUMNS * sizeof(char*));
        for (int j = 0; j < N_COLUMNS; j++) {
            retval.buffer[i][j] = (char*)malloc(COLUMN_WIDTH * sizeof(char));
        }
    }

    strcpy(retval.buffer[0][0], nlines);
    int lineI = 1;

    // Read in file
    while (fgets(line, N_COLUMNS * COLUMN_WIDTH, fp)) {
        int field = 0;
        // NOTE strtok clobbers tmp
        char* tmp = strdup(line);
        const char* tok;
        for (tok = strtok(tmp, ";");
             tok && *tok;
             tok = strtok(NULL, ";\n")) {
            strcpy(retval.buffer[lineI][field], tok);
            field++;
        }
        free(tmp);
        lineI++;
    }

    fclose(fp);

    return retval;
}