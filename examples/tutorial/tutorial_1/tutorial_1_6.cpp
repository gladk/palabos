/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <iostream>
#include <iomanip>

/* Code 1.6 in the Palabos tutorial
 */

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR plb::descriptors::D3Q19Descriptor

/// Velocity on the parabolic Poiseuille profile
T poiseuilleVelocity(plint iY, plint iZ, IncomprFlowParam<T> const& parameters) {
    T y = (T)iY / parameters.getResolution();
    T z = (T)iZ / parameters.getResolution();
    return 4.*parameters.getLatticeU() * (y-y*y) * (z-z*z);
}

/// A functional, used to initialize the velocity for the boundary conditions
template<typename T>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const {
        u[0] = poiseuilleVelocity(iY, iZ, parameters);
        u[1] = T();
        u[2] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

void channelSetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                   IncomprFlowParam<T> const& parameters,
                   OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    // Create Velocity boundary conditions
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<T>(parameters) );

    lattice.initialize();
}

void writeGifs(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter)
{
    const plint imSize = 600;
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();
    Box3D slice(0, nx-1, 0, ny-1, nz/2, nz/2);
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                               *computeVelocityNorm(lattice, slice),
                               imSize, imSize );
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    // Use the class IncomprFlowParam to convert from
    //   dimensionless variables to lattice units, in the
    //   context of incompressible flows.
    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // Reference velocity (the maximum velocity
                       //   in the Poiseuille profile) in lattice units.
            (T) 100.,  // Reynolds number
            30,        // Resolution of the reference length (channel height).
            3.,        // Channel length in dimensionless variables
            1.,        // Channel height in dimensionless variables
            1.         // Channel depth in dimensionless variables
    );
    const T imSave   = (T)0.02; // Time intervals at which to save GIF
                                //   images, in dimensionless time units.
    const T maxT     = (T)2.5;  // Total simulation time, in dimensionless
                                //   time units.

    writeLogFile(parameters, "3D Poiseuille flow");

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
             parameters.getNx(), parameters.getNy(), parameters.getNz(),
             new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>*
        //boundaryCondition = createInterpBoundaryCondition3D<T,DESCRIPTOR>();
        boundaryCondition = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    channelSetup(lattice, parameters, *boundaryCondition);

    // Main loop over time iterations.
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {
        if (iT%parameters.nStep(imSave)==0) {
            pcout << "Saving Gif at time step " << iT << endl;
            writeGifs(lattice, iT);
        }
        // Execute lattice Boltzmann iteration.
        lattice.collideAndStream();
    }

    delete boundaryCondition;
}
