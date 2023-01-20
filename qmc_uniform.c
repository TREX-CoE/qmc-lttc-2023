#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stddef.h>  // for size_t
#include "hydrogen.h"
#include "qmc_stats.h"   // for ave_error

void uniform_montecarlo(double a, size_t nmax, double *energy) {
    size_t istep;
    double norm, r[3], w;

    *energy = 0.0;
    norm = 0.0;

    for (istep = 0; istep < nmax; istep++) {
        for (int i = 0; i < 3; i++) {
            r[i] = drand48();
        }

        r[0] = -5.0 + 10.0 * r[0];
        r[1] = -5.0 + 10.0 * r[1];
        r[2] = -5.0 + 10.0 * r[2];
        w = psi(a, r);
        w = w*w;
        *energy += w * e_loc(a, r);
        norm += w;
    }
    *energy = *energy / norm;
}

int main(void) {

#define a     1.2
#define nmax  100000
#define nruns 30

    double X[nruns];
    double ave, err;

    srand48(time(NULL));

    for (size_t irun = 0; irun < nruns; irun++) {
        uniform_montecarlo(a, nmax, &X[irun]);
    }
    ave_error(X, nruns, &ave, &err);

    printf("E = %f +/- %f\n", ave, err);

    return 0;
}
