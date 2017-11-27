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

void cavitySetup( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                  IncomprFlowParam<T> const& parameters,
                  OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();

    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    setBoundaryVelocity(lattice, lattice.getBoundingBox(), Array<T,2>((T)0.,(T)0.) );
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T)1., Array<T,2>((T)0.,(T)0.) );

    T u = parameters.getLatticeU();
    setBoundaryVelocity(lattice, Box2D(1, nx-2, ny-1, ny-1), Array<T,2>(u,(T)0.) );
    initializeAtEquilibrium(lattice, Box2D(1, nx-2, ny-1, ny-1), (T)1., Array<T,2>(u,(T)0.) );

    lattice.initialize();
}

template<class BlockLatticeT>
void writeGifs(BlockLatticeT& lattice, plint iter)
{
    const plint imSize = 600;

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("uz", iter, 6),
                               *computeKineticEnergy(lattice),
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
    const plint logIt  =   100;
    const plint imSave =  1000;
    const plint maxIt  = 10000;

    writeLogFile(parameters, "2D cavity");

    // Attention: This program works only in serial, because in parallel, the two lattices
    //   refer to a different number of total threads (3*4=12 for lattice1, and 1*2=2 for lattice2).

    plint envelopeWidth = 1;
    MultiBlockLattice2D<T, DESCRIPTOR> lattice1 (
        MultiBlockManagement2D( createRegularDistribution2D (
                                  parameters.getNx(),parameters.getNy(), 3, 4),
                                  defaultMultiBlockPolicy2D().getThreadAttribution(), envelopeWidth ),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiCellAccess<T,DESCRIPTOR>(),
        new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega())
    );
    MultiBlockLattice2D<T, DESCRIPTOR> lattice2 (
        MultiBlockManagement2D( createRegularDistribution2D (
                                  parameters.getNx(),parameters.getNy(), 1, 2),
                                  defaultMultiBlockPolicy2D().getThreadAttribution(), envelopeWidth ),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiCellAccess<T,DESCRIPTOR>(),
        new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega())
    );


    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createInterpBoundaryCondition2D<T,DESCRIPTOR>();
        //boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    cavitySetup(lattice1, parameters, *boundaryCondition);
    cavitySetup(lattice2, parameters, *boundaryCondition);

    // Main loop over time iterations.
    for (plint iT=0; iT<maxIt; ++iT) {
        if (iT%logIt==0) {
            T storedAverageDensity1 = getStoredAverageDensity<T>(lattice1);
            T storedAverageDensity2 = getStoredAverageDensity<T>(lattice2);
            T storedAverageEnergy1 = getStoredAverageEnergy<T>(lattice1);
            T storedAverageEnergy2 = getStoredAverageEnergy<T>(lattice2);
            pcout << "step " << iT
                  << "; lattice 1"
                  << "; av energy="
                  << setprecision(10) << storedAverageEnergy1
                  << "; av rho="
                  << storedAverageDensity1 << endl;

            pcout << "step " << iT
                  << "; lattice 2"
                  << "; av energy="
                  << setprecision(10) << storedAverageEnergy2
                  << "; av rho="
                  << storedAverageDensity2 << endl;
        }

        if (iT%imSave==0) {
            pcout << "Saving Gif ..." << endl;
            writeGifs(lattice1, iT);
        }

        // Lattice Boltzmann iteration step.
        lattice1.collideAndStream();
        lattice2.collideAndStream();
    }

    delete boundaryCondition;
}
