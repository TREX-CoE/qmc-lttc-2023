#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stddef.h>  // for size_t
#include "hydrogen.h"
#include "qmc_stats.h"   // for ave_error

void uniform_montecarlo(double a, size_t nmax, double *energy) {
    size_t istep;
    double normalization, r[3], w, f;

    *energy = 0.0;
    normalization = 0.0;

    for (istep = 0; istep < nmax; istep++) {
        for (int i = 0; i < 3; i++) {
            r[i] = drand48();
        }

        r[0] = -5.0 + 10.0 * r[0];
        r[1] = -5.0 + 10.0 * r[1];
        r[2] = -5.0 + 10.0 * r[2];
        f = psi(a, r);
        w = f*f;
        *energy += w * e_loc(a, r);
        normalization += w;
    }
    *energy = *energy / normalization;
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
