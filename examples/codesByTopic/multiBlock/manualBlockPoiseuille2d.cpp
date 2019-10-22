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
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

/// Velocity on the parabolic Poiseuille profile
T poiseuilleVelocity(plint iY, IncomprFlowParam<T> const& parameters) {
    T y = (T)iY / parameters.getResolution();
    return 4.*parameters.getLatticeU() * (y-y*y);
}

/// Linearly decreasing pressure profile
T poiseuillePressure(plint iX, IncomprFlowParam<T> const& parameters) {
    T Lx = parameters.getNx()-1;
    T Ly = parameters.getNy()-1;
    return 8.*parameters.getLatticeNu()*parameters.getLatticeU() / (Ly*Ly) * (Lx/(T)2-(T)iX);
}

/// Convert pressure to density according to ideal gas law
T poiseuilleDensity(plint iX, IncomprFlowParam<T> const& parameters) {
    return poiseuillePressure(iX,parameters)*DESCRIPTOR<T>::invCs2 + (T)1;
}

/// A functional, used to initialize the velocity for the boundary conditions
template<typename T>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, Array<T,2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

/// A functional, used to create an initial condition for the density and velocity
template<typename T>
class PoiseuilleVelocityAndDensity {
public:
    PoiseuilleVelocityAndDensity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, T& rho, Array<T,2>& u) const {
        rho = poiseuilleDensity(iX,parameters);
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

void definePoiseuilleGeometry( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                               IncomprFlowParam<T> const& parameters,
                               OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )
{
    // Create Velocity boundary conditions
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<T>(parameters) );
    initializeAtEquilibrium (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocityAndDensity<T>(parameters) );

    lattice.initialize();
}

void writeGifs(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
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

    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 100.,  // Re
            128,        // N
            1.,        // lx
            1.         // ly 
    );

    writeLogFile(parameters, "Poiseuille flow");

    plint envelopeWidth = 1;
    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
        MultiBlockManagement2D( createRegularDistribution2D(parameters.getNx(),parameters.getNy(), 1, 1),
                                defaultMultiBlockPolicy2D().getThreadAttribution(), envelopeWidth ),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiCellAccess<T,DESCRIPTOR>(),
        new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega())
    );

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createInterpBoundaryCondition2D<T,DESCRIPTOR>();
        //boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    definePoiseuilleGeometry(lattice, parameters, *boundaryCondition);

    // Main loop over time iterations.
    for (plint iT=0; iT<1000; ++iT) {
        if (iT%100==0) {
            pcout << "step " << iT
                  << "; lattice time=" << lattice.getTimeCounter().getTime()
                  << "; t=" << iT*parameters.getDeltaT()
                  << "; av energy="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho="
                  << getStoredAverageDensity<T>(lattice) << endl;
        }

        if (iT%200==0) {
            pcout << "Saving Gif ..." << endl;
            writeGifs(lattice, iT);
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();

        Array<T,DESCRIPTOR<T>::q> populations;
        lattice.get(lattice.getNx()/2,lattice.getNy()-1).getPopulations(populations);
        pcout << populations[8] << endl;
    }

    delete boundaryCondition;
}
