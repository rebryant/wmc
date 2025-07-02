#include <stdlib.h>
#include <stdio.h>

#include "er_double.h"

void er_check(double a, double b);


int main(int argc, char *argv[]) {
    if (argc != 3) {
	fprintf(stderr, "Usage: %s A B\n", argv[0]);
	return 0;
    }
    double a = atof(argv[1]);
    double b = atof(argv[2]);
    er_check(a, b);
    return 0;
}

