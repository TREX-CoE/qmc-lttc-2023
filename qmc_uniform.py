#!/usr/bin/env python3

from hydrogen  import *
from qmc_stats import *
import numpy as np

def MonteCarlo(a, nmax):
     energy = 0.
     normalization = 0.

     for istep in range(nmax):
          #r = np.random.uniform(-5., 5., (3))
          R = 5.
          phi = np.random.rand()*2*np.pi
          costheta = np.random.rand()*2 - 1.0
          u = np.random.rand()

          theta = np.arccos( costheta )
          r = R * np.cbrt( u )

          # Spherical distribution
          x = r * np.sin( theta) * np.cos( phi )
          y = r * np.sin( theta) * np.sin( phi )
          z = r * np.cos( theta )
          r = np.array([x,y,z],dtype=np.Float64)


          w = psi(a,r)
          w = w*w

          energy        += w * e_loc(a,r)
          normalization += w

     return energy / normalization

a    = 1.2
nmax = 100000

X = [MonteCarlo(a,nmax) for i in range(30)]
E, deltaE = ave_error(X)

print(f"E = {E} +/- {deltaE}")
