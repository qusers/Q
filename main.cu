#include "system.h"

#include <math.h>
#include <stdio.h>

int main(int argc, char *argv[]){
    srand(1);

    if (argc < 2) {
        printf(">>> FATAL: Input file folder expected. Exiting...\n");
        exit(EXIT_FAILURE);
    }

    strcpy(base_folder, argv[1]);

    init_variables();

    for (int i = 0; i < 11; i++) {
        calc_integration(i);
    }

    clean_variables();

    exit(EXIT_SUCCESS);
}