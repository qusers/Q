#include "system.h"

#include <math.h>
#include <stdio.h>

int main(){
    srand(1);

    init_variables();

    for (int i = 0; i < 10001; i++) {
        calc_integration(i);
    }

    clean_variables();

    exit(EXIT_SUCCESS);
}