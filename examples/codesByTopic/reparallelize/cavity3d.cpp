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
    initializeAtEquilibrium(lattice, everythingButTopLid, (T) 1., Array<T,3>((T)0.,u/(T)3.,(T)0.) );
    initializeAtEquilibrium(lattice, topLid, (T) 1., Array<T,3>(u,(T)0.,u) );
    setBoundaryVelocity(lattice, topLid, Array<T,3>(u,(T)0.,u) );

    lattice.initialize();
}

template<class BlockLatticeT>
void writeGifs(BlockLatticeT& lattice,
               IncomprFlowParam<T> const& parameters, plint iter)
{
    const plint imSize = 600;
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();
    //const plint zComponent = 2;

    Box3D slice(0, nx-1, 0, ny-1, nz/2, nz/2);
    ImageWriter<T> imageWriter("leeloo");

    auto_ptr<MultiScalarField3D<T> > velNorm(computeVelocityNorm(lattice));
    auto_ptr<MultiScalarField3D<T> > velNorm2 = redistribute (
            *velNorm, createRegularDistribution3D (
                        lattice.getBoundingBox(),
                        // Re-parallelize in z-direction.
                        1, 1, global::mpi().getSize() ) );

    /*
    imageWriter.writeScaledGif( createFileName("uz", iter, 6),
                                *computeVelocityComponent (lattice, slice, zComponent),
                                imSize, imSize );
                                */

    imageWriter.writeScaledGif( createFileName("uNorm", iter, 6),
                                *extractSubDomain(*velNorm2, slice),
                                //*computeVelocityNorm (lattice, slice),
                                imSize, imSize );
    /*
    imageWriter.writeScaledGif( createFileName("omega", iter, 6),
                                *computeNorm(*computeVorticity (
                                        *computeVelocity(lattice) ), slice ),
                                imSize, imSize );
                                */
}

template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
    vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
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
    const T logT     = (T)1/(T)100;
    //const T imSave   = (T)1/(T)40;
    const T vtkSave  = (T)1;
    const T maxT     = (T)10.1;

    writeLogFile(parameters, "3D diagonal cavity");

    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();

    std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > originLattice (
        new MultiBlockLattice3D<T,DESCRIPTOR> (
            MultiBlockManagement3D (
                // Parallelize in x-direction
                createRegularDistribution3D (
                    parameters.getNx(), parameters.getNy(), parameters.getNz(),
                    global::mpi().getSize(), 1, 1),              // envelopeWidth
                defaultMultiBlockPolicy3D().getThreadAttribution(), 1 ),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
            new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) )
    );


    //lattice->setInternalTypeOfModification(modif::dataStructure);

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createInterpBoundaryCondition3D<T,DESCRIPTOR>();
        //= createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    cavitySetup(*originLattice, parameters, *boundaryCondition);
    Box3D centerBox(nx/3, 2*nx/3, ny/3, 2*ny/3, nz/3, 2*nz/3);
    auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > intermediateLattice
     /*
        = redistribute (
                *originLattice,
                createRegularDistribution3D (
                    originLattice->getBoundingBox(),
                    // Re-parallelize in y-direction.
                    1, global::mpi().getSize(), 1 )
          );
          */
        = except(*originLattice, centerBox);
    //auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > lattice
    //    = extend(*intermediateLattice, Box3D(nx/3+nx/12, 2*nx/3-nx/12, 0, ny-1, 0, nz-1));

    defineDynamics(*intermediateLattice, centerBox.enlarge(1), new BounceBack<T,DESCRIPTOR>(1.));
    //auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > lattice = reparallelize(*intermediateLattice);
    MultiBlockLattice3D<T,DESCRIPTOR> dummyLattice(nx/2, ny/2, nz/2, new BounceBack<T,DESCRIPTOR>(1.));
    auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > lattice = align(*intermediateLattice, dummyLattice);

    T previousIterationTime = T();
    // Loop over main time iteration.
    for (plint iT=0; iT<parameters.nStep(maxT); ++iT) {
        global::timer("mainLoop").restart();

        if (iT%100==0) {
        //if (iT%parameters.nStep(imSave)==0) {
            pcout << "Writing Gif ..." << endl;
            writeGifs(*lattice, parameters, iT);
        }

        if (iT%parameters.nStep(vtkSave)==0 && iT>0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(*lattice, parameters, iT);
        }

        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; t=" << iT*parameters.getDeltaT();
        }

        // Execute a time iteration.
        lattice->collideAndStream();

        // Access averages from internal statistics ( their value is defined
        //   only after the call to lattice.collideAndStream() )
        if (iT%parameters.nStep(logT)==0) {
            pcout << "; av energy="
                  << setprecision(10) << getStoredAverageEnergy<T>(*lattice)
                  << "; av rho="
                  << setprecision(10) << getStoredAverageDensity<T>(*lattice) << endl;
            pcout << "Time spent during previous iteration: "
                  << previousIterationTime << endl;
        }

        previousIterationTime = global::timer("mainLoop").stop();
    }

    delete boundaryCondition;
}
