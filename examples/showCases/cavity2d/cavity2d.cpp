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

/** \file
  * Flow in a lid-driven 2D cavity. The cavity is square and has no-slip walls,
  * except for the top wall which is driven to the right with a constant
  * velocity. The benchmark is challenging because of the velocity
  * discontinuities on corner nodes. The code on the other hand is very simple.
  * It could for example be used as a first example, to get familiar with Palabos.
  **/

#include "palabos2D.h"
#ifndef PLB_PRECOMPILED // Unless precompiled version is used,
  #include "palabos2D.hh"   // include full template code
#endif
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
void writeGif(BlockLatticeT& lattice, plint iter)
{
    const plint imSize = 600;

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("uNorm", iter, 6),
                               *computeVelocityNorm(lattice),
                               imSize, imSize );
    imageWriter.writeScaledGif(createFileName("logUnorm", iter, 6),
                               *computeLog(*add((T)1.e-8,*computeVelocityNorm(lattice))),
                               imSize, imSize );
}

template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity", dx/dt);
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 100.,  // Re
            128,       // N
            1.,        // lx
            1.         // ly
    );
    const T logT     = (T)0.1;
    const T imSave   = (T)0.2;
    const T vtkSave  = (T)1.;
    const T maxT     = (T)10.1;

    writeLogFile(parameters, "2D cavity");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
              parameters.getNx(), parameters.getNy(),
              new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        //boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();
        boundaryCondition = createInterpBoundaryCondition2D<T,DESCRIPTOR>();

    cavitySetup(lattice, parameters, *boundaryCondition);

    T previousIterationTime = T();

    // Main loop over time iterations.
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {
        global::timer("mainLoop").restart();

        if (iT%parameters.nStep(imSave)==0 && iT>0) {
            pcout << "Saving Gif ..." << endl;
            writeGif(lattice, iT);
            pcout << endl;
        }

        if (iT%parameters.nStep(vtkSave)==0 && iT>0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(lattice, parameters, iT);
        }

        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; t=" << iT*parameters.getDeltaT();
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();

        if (iT%parameters.nStep(logT)==0) {
            pcout << "; av energy="
                  << setprecision(10) << getStoredAverageEnergy(lattice)
                  << "; av rho="
                  << getStoredAverageDensity(lattice) << endl;
            pcout << "Time spent during previous iteration: "
                  << previousIterationTime << endl;
        }

        previousIterationTime = global::timer("mainLoop").stop();
    }

    delete boundaryCondition;
}
