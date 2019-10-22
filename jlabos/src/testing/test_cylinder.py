#!/usr/bin/python

#  This file is part of the Palabos library.
#  Copyright (C) 2011-2017 FlowKit Sarl
#  Route d'Oron 2
#  1010 Lausanne, Switzerland
#  E-mail contact: contact@flowkit.com
#
#  The library Palabos is free software: you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.

import sys; sys.path.append('..')
from pypal import *
import numpy

def significant(val):
    return round(val,10)

def cylinderFlow2D(nx, ny, numIter, floatT):
    """
    Simulates a flow in a lid-driven 2D cavity, and returns the average energy
    and density. This function was written as a regression test for the
    Palabos library.

    Different size, but square shape:
    >>> [cylinderFlow2D(n,n, 100, float) for n in range(60,140,19)]
    [\
['Energy OK', True, 'Density OK', True, 'Av. energy', '4.324600e-06', 'Av. density', '9.999913e-01'], \
['Energy OK', True, 'Density OK', True, 'Av. energy', '3.504000e-06', 'Av. density', '9.999951e-01'], \
['Energy OK', True, 'Density OK', True, 'Av. energy', '2.989000e-06', 'Av. density', '9.999969e-01'], \
['Energy OK', True, 'Density OK', True, 'Av. energy', '2.635400e-06', 'Av. density', '9.999978e-01'], \
['Energy OK', True, 'Density OK', True, 'Av. energy', '2.375300e-06', 'Av. density', '9.999984e-01']]

    Non-square aspect ratio:
    >>> cylinderFlow2D(100,200, 200, float)
    ['Energy OK', True, 'Density OK', True, 'Av. energy', '2.279700e-06', 'Av. density', '9.999971e-01']

    Single-precision floating point:
    >>> cylinderFlow2D(100,200, 200, numpy.float32)
    ['Energy OK', True, 'Density OK', True, 'Av. energy', '2.279700e-06', 'Av. density', '9.999971e-01']
    """

    uMax     = 0.02          # Maximum velocity in lattice units.
    Re       = 100.          # Reynolds number (wrt ny).
    nu       = uMax*ny / Re  # Viscosity in lattice units.
    tau      = 3*nu+0.5      # Relaxation time.

    # Allocate memory for the block-lattice, and assign the dynamics.
    lattice = Block2D(nx, ny, D2Q9, BGK(omega=floatT(1./tau)))

    # Use a regularized boundary condition.
    boundary.regularized().defineVelocityBC(lattice)
    # Default all external walls to zero velocity.
    lattice.setBoundaryVelocity([0.,0.])
    # Top lid has constant in-plane velocity.
    lattice[1:nx-2,ny-1].setBoundaryVelocity([uMax,0.])

    lattice.initializeAtEquilibrium(1., [0.,0.])

    for i in range(0,numIter):
        # This is the real computation: evolving through a collision-streaming cycle.
        lattice.collideAndStream()
    
    averageEnergy1 = significant( lattice.averageEnergy() )
    averageEnergy2 = significant( lattice.kineticEnergy().average() )

    averageDensity1 = significant( lattice.averageDensity() )
    averageDensity2 = significant( lattice.density().average() )

    return [ 'Energy OK', averageEnergy1 == averageEnergy2,
             'Density OK', averageDensity1 == averageDensity2,
             'Av. energy', '%.6e'%averageEnergy1,
             'Av. density', '%.6e'%averageDensity1 ]


def cylinderFlow3D(nx, ny, nz, numIter, floatT):
    """
    Simulates a flow in a lid-driven 3D cavity, and returns the average energy
    and density. This function was written as a regression test for the
    Palabos library.

    Different size, but square shape:
    >>> [cylinderFlow3D(n,n,n, 30, float) for n in range(20,60,19)]
    [\
['Energy OK', True, 'Density OK', True, 'Av. energy', '8.491700e-06', 'Av. density', '9.999741e-01'], \
['Energy OK', True, 'Density OK', True, 'Av. energy', '4.915800e-06', 'Av. density', '9.999932e-01'], \
['Energy OK', True, 'Density OK', True, 'Av. energy', '3.531500e-06', 'Av. density', '9.999969e-01']]

    Non-square aspect ratio:
    >>> cylinderFlow3D(20,31,42, 30, float)
    ['Energy OK', True, 'Density OK', True, 'Av. energy', '4.277200e-06', 'Av. density', '9.999879e-01']

    Single-precision floating point:
    >>> cylinderFlow3D(20,31,42, 30, numpy.float32)
    ['Energy OK', True, 'Density OK', True, 'Av. energy', '4.277200e-06', 'Av. density', '9.999879e-01']
    """

    uMax     = 0.02          # Maximum velocity in lattice units.
    Re       = 100.          # Reynolds number (wrt ny).
    nu       = uMax*ny / Re  # Viscosity in lattice units.
    tau      = 3*nu+0.5      # Relaxation time.

    # Allocate memory for the block-lattice, and assign the dynamics.
    lattice = Block3D(nx, ny, nz, D3Q19, BGK(omega=floatT(1./tau)))

    # Use a regularized boundary condition.
    boundary.regularized().defineVelocityBC(lattice)
    # Default all external walls to zero velocity.
    lattice.setBoundaryVelocity([0.,0.,0.])
    # Top lid has constant in-plane velocity.
    lattice[1:nx-2,1:ny-2,nz-1].setBoundaryVelocity([uMax,0.,0.])

    lattice.initializeAtEquilibrium(1., [0.,0.,0.])

    for i in range(0,numIter):
        # This is the real computation: evolving through a collision-streaming cycle.
        lattice.collideAndStream()
    
    averageEnergy1 = significant( lattice.averageEnergy() )
    averageEnergy2 = significant( lattice.kineticEnergy().average() )

    averageDensity1 = significant( lattice.averageDensity() )
    averageDensity2 = significant( lattice.density().average() )

    return [ 'Energy OK', averageEnergy1 == averageEnergy2,
             'Density OK', averageDensity1 == averageDensity2,
             'Av. energy', '%.6e'%averageEnergy1,
             'Av. density', '%.6e'%averageDensity1 ]


def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()

