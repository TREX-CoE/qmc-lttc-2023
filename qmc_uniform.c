#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "hydrogen.h"
#include "qmc_stats.h"   // for ave_error

void uniform_montecarlo(double a, long long int nmax, double *energy) {
    long long int istep;
    double norm, r[3], w;

    *energy = 0.0;
    norm = 0.0;

    for (istep = 0; istep < nmax; istep++) {
        for (int i = 0; i < 3; i++) {
            r[i] = (double)rand() / RAND_MAX;
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
    double a = 1.2;
    long long int nmax = 100000;
    int nruns = 30;

    double X[nruns];
    double ave, err;

    for (int irun = 0; irun < nruns; irun++) {
        uniform_montecarlo(a, nmax, &X[irun]);
    }
    ave_error(X, nruns, &ave, &err);

    printf("E = %f +/- %f\n", ave, err);

    return 0;
}
