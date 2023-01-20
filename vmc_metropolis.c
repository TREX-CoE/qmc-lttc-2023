#include <stdio.h>
#include <stdlib.h>
#include <stddef.h> // for size_t
#include <math.h>
#include <time.h>
#include "hydrogen.h"
#include "qmc_stats.h"

void variational_montecarlo(double a, size_t nmax, double dt, 
                            double *energy, double *accep)
{
    double chi[3], d2_old, prod, u;
    double psi_old, psi_new, d2_new, argexpo, q;
    double r_old[3], r_new[3];
    double d_old[3], d_new[3];
    size_t istep, n_accep = 0;

    double sq_dt = sqrt(dt);

    // Initialization
    *energy = 0.0;

    random_gauss(r_old, 3);

    drift(a, r_old, d_old);
    d2_old = d_old[0]*d_old[0] + d_old[1]*d_old[1] + d_old[2]*d_old[2];

    psi_old = psi(a, r_old);

    for (istep = 0; istep < nmax; istep++) {
        *energy += e_loc(a, r_old);

        random_gauss(chi, 3);
        for (int i = 0; i < 3; i++) {
            r_new[i] = r_old[i] + dt*d_old[i] + chi[i]*sq_dt;
        }

        drift(a, r_new, d_new);
        d2_new = d_new[0]*d_new[0] + d_new[1]*d_new[1] + d_new[2]*d_new[2];

        psi_new = psi(a, r_new);

        // Metropolis
        prod = (d_new[0] + d_old[0])*(r_new[0] - r_old[0]) +
               (d_new[1] + d_old[1])*(r_new[1] - r_old[1]) +
               (d_new[2] + d_old[2])*(r_new[2] - r_old[2]);
        argexpo = 0.5 * (d2_new - d2_old)*dt + prod;

        q = psi_new / psi_old;
        q = exp(-argexpo) * q*q;

        u = drand48();

        if (u <= q) {
            n_accep++;
            for (int i = 0; i < 3; i++) {
                r_old[i] = r_new[i];
                d_old[i] = d_new[i];
            }
            d2_old = d2_new;
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
        variational_montecarlo(a, nmax, dt, &X[irun], &Y[irun]);
    }

    ave_error(X, nruns, &ave, &err);
    printf("E = %f +/- %f\n", ave, err);

    ave_error(Y, nruns, &ave, &err);
    printf("A = %f +/- %f\n", ave, err);

    return 0;
}
