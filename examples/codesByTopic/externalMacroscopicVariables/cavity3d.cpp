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
  * Flow in a lid-driven 3D cavity. The cavity is square and has no-slip walls,
  * except for the top wall which is diagonally driven with a constant
  * velocity. The benchmark is challenging because of the velocity
  * discontinuities on corner nodes.
  **/

#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

void cavitySetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                  IncomprFlowParam<T> const& parameters,
                  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    Box3D topLid = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);
    Box3D everythingButTopLid = Box3D(0, nx-1, 0, ny-2, 0, nz-1);

    // All walls implement a Dirichlet velocity condition.
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    T u = std::sqrt((T)2)/(T)2 * parameters.getLatticeU();
    initializeAtEquilibrium(lattice, everythingButTopLid, (T) 1., Array<T,3>((T)0.,(T)0.,(T)0.) );
    initializeAtEquilibrium(lattice, topLid, (T) 1., Array<T,3>(u,(T)0.,u) );
    setBoundaryVelocity(lattice, topLid, Array<T,3>(u,(T)0.,u) );

    lattice.initialize();
}

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    //defaultMultiBlockPolicy3D().setNumProcesses(8);

    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 100.,  // Re
            100,        // N
            1.,        // lx
            1.,        // ly
            1.         // lz
    );
    const T maxT     = (T)10.1;

    writeLogFile(parameters, "3D diagonal cavity");

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(), parameters.getNz(),
            new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    MultiBlockLattice3D<T, DESCRIPTOR> lattice2 (
            parameters.getNx(), parameters.getNy(), parameters.getNz(),
            new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        //= createInterpBoundaryCondition3D<T,DESCRIPTOR>();
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    cavitySetup(lattice, parameters, *boundaryCondition);
    cavitySetup(lattice2, parameters, *boundaryCondition);

    MultiScalarField3D<T> rhoBar(lattice2);
    MultiTensorField3D<T,3> j(lattice2);
    std::vector<MultiBlock3D*> lattice2RhoBarJparam;
    lattice2RhoBarJparam.push_back(&lattice2);
    lattice2RhoBarJparam.push_back(&rhoBar);
    lattice2RhoBarJparam.push_back(&j);
    computeRhoBarJ(lattice2, rhoBar, j, lattice2.getBoundingBox());

    integrateProcessingFunctional( new ExternalRhoJcollideAndStream3D<T,DESCRIPTOR>,
                                   lattice2.getBoundingBox(), lattice2RhoBarJparam, -1 );

    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    Box3D slice(0.5*nx-1, 0.5*nx+1, 0, ny-1, 0, nz-1);

    // Loop over main time iteration.
    global::timer("iterations").start();
    for (plint iT=0; iT<parameters.nStep(maxT); ++iT) {

        // Execute a time iteration.
        lattice.collideAndStream();
        lattice2.executeInternalProcessors(-1);
        lattice2.executeInternalProcessors();
        computeRhoBarJ(lattice2, rhoBar, j, lattice2.getBoundingBox());

        if (iT%100==0) {
            pcout << "Energy 1: " << computeAverageEnergy(lattice) << std::endl;
            pcout << "Energy 2: " << computeAverageEnergy(lattice2) << std::endl;
            pcout << global::timer("iterations").getTime()/(iT+1) << std::endl;

            T dx = parameters.getDeltaX();
            T dt = parameters.getDeltaT();
            VtkImageOutput3D<T> vtkOut(createFileName("diff_", iT, 6), dx);
            vtkOut.writeData<float>(
                    *subtract(*computeVelocityNorm(lattice, slice),*computeVelocityNorm(lattice2, slice)), 
                    "velocityNormDiff", dx/dt);
        }

    }

    delete boundaryCondition;
}
