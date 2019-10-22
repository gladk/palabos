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
#include <sstream>
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


void getParametersFromCommandLine(int argc, char* argv[], T& Re, plint& N) {
    // Make sure the command line has two paramters.
    if (argc != 3) {
        pcout << "Error; Correct syntax is \" " << argv[0] << " Re N\"" << endl;
        exit(-1);
    } 

    // Read Reynolds number and resolution from command line.
    stringstream ReStr, NStr;
    ReStr << argv[1]; ReStr >> Re;
    NStr  << argv[2]; NStr >> N;
}

void getParametersFromParamFile(T& Re, plint& N) {
    plb_ifstream ifile("parameters.dat");
    ifile >> Re >> N;
    global::mpi().bCast(&Re, 1);
    global::mpi().bCast(&N, 1);
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    T    Re; // Reynolds number.
    plint N;  // Resolution.

    //getParametersFromCommandLine(argc, argv, Re, N);
    getParametersFromParamFile(Re, N);

    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 100.,  // Re
            128,        // N
            1.,        // lx
            1.         // ly 
    );
    const T logT     = (T)0.1;
    const T imSave   = (T)0.2;
    const T maxT     = (T)10.1;

    writeLogFile(parameters, "2D cavity");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
              parameters.getNx(), parameters.getNy(),
              new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        //boundaryCondition = createInterpBoundaryCondition2D<T,DESCRIPTOR>();
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    cavitySetup(lattice, parameters, *boundaryCondition);

    // Main loop over time iterations.
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {
        if ((iT+1)%parameters.nStep(logT)==0) {
            pcout << computeAverageDensity(lattice) << endl;
            pcout << computeAverageEnergy(lattice) << endl;
        }
        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; lattice time=" << lattice.getTimeCounter().getTime()
                  << "; t=" << iT*parameters.getDeltaT()
                  << "; av energy="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                   << "; av rho="
                 << getStoredAverageDensity<T>(lattice) << endl;
        }

        if (iT%parameters.nStep(imSave)==0 && iT>0) {
            pcout << "Saving Gif ..." << endl;
            writeGifs(lattice, iT);
            pcout << endl;
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();
    }

    delete boundaryCondition;
}
