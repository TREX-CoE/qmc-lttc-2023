#include <stdio.h>
#include <math.h>
#include "hydrogen.h"

#define NPOINTS  50
#define NEXPO     6

int main() {

    double x[NPOINTS], energy, dx, r[3];
    double a[NEXPO] = { 0.1, 0.2, 0.5, 1.0, 1.5, 2.0 };
    int i, j;

    dx = 10.0/(NPOINTS-1);
    for (i = 0; i < NPOINTS; i++) {
        x[i] = -5.0 + i*dx;
    }

    // TODO
    return 0;
}

#include <stdio.h>
#include <math.h>
#include "hydrogen.h"

#define NPOINTS  50
#define NEXPO     6

int main() {

    double x[NPOINTS], energy, dx, r[3];
    double a[NEXPO] = { 0.1, 0.2, 0.5, 1.0, 1.5, 2.0 };
    int i, j;

    dx = 10.0/(NPOINTS-1);
    for (i = 0; i < NPOINTS; i++) {
        x[i] = -5.0 + i*dx;
    }

    for (i = 0; i < 3; i++) {
        r[i] = 0.0;
    }

    for (j = 0; j < NEXPO; j++) {
        printf("# a=%f\n", a[j]);
        for (i = 0; i < NPOINTS; i++) {
            r[0] = x[i];
            energy = e_loc(a[j], r);
            printf("%f %f\n", x[i], energy);
        }
        printf("\n\n");
    }
    return 0;
}
