#include <stdio.h>
#include <math.h>
#include "hydrogen.h"

#define NPOINTS  50
#define NEXPO     6

int main() {

    double x[NPOINTS], energy, dx, r[3], delta, norm, w;
    double a[NEXPO] = { 0.1, 0.2, 0.5, 1.0, 1.5, 2.0 };
    double energy2, e_tmp, s2;

    dx = 10.0/(NPOINTS-1);
    for (int i = 0; i < NPOINTS; i++) {
        x[i] = -5.0 + i*dx;
    }

    delta = dx*dx*dx;
    for (int i = 0; i < 3; i++) {
        r[i] = 0.0;
    }

    for (int j = 0; j < NEXPO; j++) {
        energy  = 0.0;
        energy2 = 0.0;
        norm    = 0.0;

        for (int i = 0; i < NPOINTS; i++) {
            r[0] = x[i];

            for (int k = 0; k < NPOINTS; k++) {
                r[1] = x[k];

                for (int l = 0; l < NPOINTS; l++) {
                    r[2] = x[l];

                    w = psi(a[j], r);
                    w = w*w*delta;

                    e_tmp = e_loc(a[j], r);

                    energy  += w * e_tmp;
                    energy2 += w * e_tmp * e_tmp;
                    norm    += w;
                }
            }
        }
        energy  = energy/norm;
        energy2 = energy2/norm;
        s2 = energy2 - energy*energy;
        printf("a = %f    E = %f    s2 = %f\n", a[j], energy, s2);
    }
}
