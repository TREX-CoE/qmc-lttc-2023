#include <stdio.h>
#include <stdlib.h>
#include <stddef.h> // for size_t
#include <math.h>
#include <time.h>
#include "hydrogen.h"
#include "qmc_stats.h"

void metropolis_montecarlo(double a, size_t nmax, double dt,
                           double *energy, double *accep)
{
    double r_old[3], r_new[3], psi_old, psi_new, v, ratio;
    size_t n_accep = 0;

    *energy = 0.0;

    for (int i = 0; i < 3; i++) {
        r_old[i] = dt * (2.0*drand48() - 1.0);
    }
    psi_old = psi(a, r_old);

    for (size_t istep = 0; istep < nmax; istep++) {
        *energy += e_loc(a, r_old);

        for (int i = 0; i < 3; i++) {
            r_new[i] = r_old[i] + dt * (2.0*drand48() - 1.0);
        }

        psi_new = psi(a, r_new);

        ratio = pow(psi_new / psi_old,2);
        v = drand48();

        if (v <= ratio) {
            n_accep++;
            for (int i = 0; i < 3; i++) {
                r_old[i] = r_new[i];
            }
            psi_old = psi_new;
        }
    }
    *energy = *energy / (double) nmax;
    *accep = (double) n_accep / (double) nmax;
}

int main(void) {

#define a      1.2
#define nmax   100000
#define dt     1.0
#define nruns  30

    double X[nruns];
    double Y[nruns];
    double ave, err;

    srand48(time(NULL));

    for (size_t irun = 0; irun < nruns; irun++) {
        metropolis_montecarlo(a, nmax, dt, &X[irun], &Y[irun]);
    }

    ave_error(X, nruns, &ave, &err);
    printf("E = %f +/- %f\n", ave, err);

    ave_error(Y, nruns, &ave, &err);
    printf("A = %f +/- %f\n", ave, err);

    return 0;
}
