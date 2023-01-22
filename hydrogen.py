#!/usr/bin/env python3
import numpy as np

def potential(r):
    distance = np.sqrt(np.dot(r,r))
    if distance == 0:
        print("potential at r=0 diverges")
        return inf
    return -1. / distance


def test_potential():
    expected_output = -1./15.
    for r in [( 2., 5., 14.), (5., 14., 2.), 
              (-2., 5.,-14.), (5.,-14.,-2.), 
              ( 0., 9.,-12.), (9.,-12., 0.)]:
        assert potential(r) == expected_output

    r = (0., 0., 0.)
    assert potential(r) == inf

    print("potential ok")

if __name__ == "__main__":
    test_potential()

def psi(a, r):
    return np.exp(-a*np.sqrt(np.dot(r,r)))

def kinetic(a,r):
    distance = np.sqrt(np.dot(r,r))
    assert (distance > 0.)

    return a * (1./distance - 0.5 * a)

def e_loc(a,r):
    return kinetic(a,r) + potential(r)

def drift(a,r):
   ar_inv = -a/np.sqrt(np.dot(r,r))
   return r * ar_inv
