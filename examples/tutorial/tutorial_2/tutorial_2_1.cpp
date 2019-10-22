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

#include "palabos2D.h"
#include "palabos2D.hh"
#include <vector>
#include <iostream>
#include <iomanip>

/* Code 2.1 in the Palabos tutorial
 */

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor

/// Velocity on the parabolic Poiseuille profile
T poiseuilleVelocity(plint iY, IncomprFlowParam<T> const& parameters) {
    T y = (T)iY / parameters.getResolution();
    return 4.*parameters.getLatticeU() * (y-y*y);
}

/// A functional, used to initialize the velocity for the boundary conditions
template<typename T>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    /// This version of the operator returns the velocity only,
    ///    to instantiate the boundary condition.
    void operator()(plint iX, plint iY, Array<T,2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
    /// This version of the operator returns also a constant value for
    ///    the density, to create the initial condition.
    void operator()(plint iX, plint iY, T& rho, Array<T,2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
        rho  = (T)1;
    }
private:
    IncomprFlowParam<T> parameters;
};

void channelSetup (
        BlockLattice2D<T,DESCRIPTOR>& lattice,
        IncomprFlowParam<T> const& parameters,
        OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )
{
    // Create Velocity boundary conditions.
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    // Specify the boundary velocity.
    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<T>(parameters) );

    // Create the initial condition.
    initializeAtEquilibrium (
         lattice, lattice.getBoundingBox(), PoiseuilleVelocity<T>(parameters) );

    lattice.initialize();
}

void writeGifs(BlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
{
    const plint imSize = 600;
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                               *computeVelocityNorm(lattice),
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
            100,       // Resolution of the reference length (channel height).
            2.,        // Channel length in dimensionless variables
            1.         // Channel height in dimensionless variables
    );
    const T imSave   = (T)0.1;  // Time intervals at which to save GIF
                                //   images, in dimensionless time units.
    const T maxT     = (T)3.1;  // Total simulation time, in dimensionless
                                //   time units.

    writeLogFile(parameters, "Poiseuille flow");

    BlockLattice2D<T, DESCRIPTOR> lattice (
             parameters.getNx(), parameters.getNy(),
             new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    channelSetup(lattice, parameters, *boundaryCondition);

    // Main loop over time iterations.
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {
        if (iT%parameters.nStep(imSave)==0 && iT>0) {
            pcout << "Saving Gif at time step " << iT << endl;
            writeGifs(lattice, iT);
        }
        // Execute lattice Boltzmann iteration.
        lattice.collideAndStream();
    }

    delete boundaryCondition;
}
