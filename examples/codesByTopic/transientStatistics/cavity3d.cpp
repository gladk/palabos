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
  *
  * This code demonstrates the use of transient statistics.
  *
  **/

#include "palabos3D.h"
#include "palabos3D.hh"
#include "palabos2D.h"
#include "palabos2D.hh"
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

    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, lattice.getBoundingBox(), boundary::dirichlet);

    T u = std::sqrt((T)2)/(T)2 * parameters.getLatticeU();
    initializeAtEquilibrium(lattice, everythingButTopLid, (T) 1., Array<T,3>((T)0.,(T)0.,(T)0.) );
    initializeAtEquilibrium(lattice, topLid, (T) 1., Array<T,3>(u,(T)0.,u) );
    setBoundaryVelocity(lattice, topLid, Array<T,3>(u,(T)0.,u) );

    lattice.initialize();
}

void initializeTransientStatistics(TransientStatistics3D<T,DESCRIPTOR>& transientStatistics)
{
    // Fields and operations must be registered in the transient statistics data structure before used.
    //
    // The allowed fields are among:
    // "velocityX", "velocityY", "velocityZ", "velocityNorm", "pressure", "vorticityX", "vorticityY", "vorticityZ", "vorticityNorm"
    //
    // The allowed operations are:
    // "min", "max", "ave" (for mean value), "rms", "dev" (for standard deviation).

    // Here we register all possible transient statistics fields and operations for demonstration
    // purposes. The user might choose a part of them in her specific application.

    (void) transientStatistics.registerFieldOperation("velocityX", "min");
    (void) transientStatistics.registerFieldOperation("velocityX", "max");
    (void) transientStatistics.registerFieldOperation("velocityX", "ave");
    (void) transientStatistics.registerFieldOperation("velocityX", "rms");
    (void) transientStatistics.registerFieldOperation("velocityX", "dev");

    (void) transientStatistics.registerFieldOperation("velocityY", "min");
    (void) transientStatistics.registerFieldOperation("velocityY", "max");
    (void) transientStatistics.registerFieldOperation("velocityY", "ave");
    (void) transientStatistics.registerFieldOperation("velocityY", "rms");
    (void) transientStatistics.registerFieldOperation("velocityY", "dev");

    (void) transientStatistics.registerFieldOperation("velocityZ", "min");
    (void) transientStatistics.registerFieldOperation("velocityZ", "max");
    (void) transientStatistics.registerFieldOperation("velocityZ", "ave");
    (void) transientStatistics.registerFieldOperation("velocityZ", "rms");
    (void) transientStatistics.registerFieldOperation("velocityZ", "dev");

    (void) transientStatistics.registerFieldOperation("velocityNorm", "min");
    (void) transientStatistics.registerFieldOperation("velocityNorm", "max");
    (void) transientStatistics.registerFieldOperation("velocityNorm", "ave");
    (void) transientStatistics.registerFieldOperation("velocityNorm", "rms");
    (void) transientStatistics.registerFieldOperation("velocityNorm", "dev");

    (void) transientStatistics.registerFieldOperation("pressure", "min");
    (void) transientStatistics.registerFieldOperation("pressure", "max");
    (void) transientStatistics.registerFieldOperation("pressure", "ave");
    (void) transientStatistics.registerFieldOperation("pressure", "rms");
    (void) transientStatistics.registerFieldOperation("pressure", "dev");

    (void) transientStatistics.registerFieldOperation("vorticityX", "min");
    (void) transientStatistics.registerFieldOperation("vorticityX", "max");
    (void) transientStatistics.registerFieldOperation("vorticityX", "ave");
    (void) transientStatistics.registerFieldOperation("vorticityX", "rms");
    (void) transientStatistics.registerFieldOperation("vorticityX", "dev");

    (void) transientStatistics.registerFieldOperation("vorticityY", "min");
    (void) transientStatistics.registerFieldOperation("vorticityY", "max");
    (void) transientStatistics.registerFieldOperation("vorticityY", "ave");
    (void) transientStatistics.registerFieldOperation("vorticityY", "rms");
    (void) transientStatistics.registerFieldOperation("vorticityY", "dev");

    (void) transientStatistics.registerFieldOperation("vorticityZ", "min");
    (void) transientStatistics.registerFieldOperation("vorticityZ", "max");
    (void) transientStatistics.registerFieldOperation("vorticityZ", "ave");
    (void) transientStatistics.registerFieldOperation("vorticityZ", "rms");
    (void) transientStatistics.registerFieldOperation("vorticityZ", "dev");

    (void) transientStatistics.registerFieldOperation("vorticityNorm", "min");
    (void) transientStatistics.registerFieldOperation("vorticityNorm", "max");
    (void) transientStatistics.registerFieldOperation("vorticityNorm", "ave");
    (void) transientStatistics.registerFieldOperation("vorticityNorm", "rms");
    (void) transientStatistics.registerFieldOperation("vorticityNorm", "dev");

    // Next, we need to initialize all transient statistics.
    // The statistics are initialized with the value of the current field
    // for the operations "min", "max", and "ave". They are initialized
    // with the absolute value of the field for "rms", and with zero for "dev".
    // The function "initialize" must be called once, and can be called
    // at any point in the simulation.
    transientStatistics.initialize();
}

template<class BlockLatticeT>
void writeGifs(BlockLatticeT& lattice,
               IncomprFlowParam<T> const& parameters, plint iter)
{
    const plint imSize = 600;
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    const plint zComponent = 2;

    Box3D slice(0, nx-1, 0, ny-1, nz/2, nz/2);
    ImageWriter<T> imageWriter("leeloo");

    imageWriter.writeScaledGif( createFileName("uz", iter, 6),
                                *computeVelocityComponent (lattice, slice, zComponent),
                                imSize, imSize );

    imageWriter.writeScaledGif( createFileName("uNorm", iter, 6),
                                *computeVelocityNorm (lattice, slice),
                                imSize, imSize );
    imageWriter.writeScaledGif( createFileName("omega", iter, 6),
                                *computeNorm(*computeVorticity (
                                        *computeVelocity(lattice) ), slice ),
                                imSize, imSize );
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

    plb_ofstream ofile("velocity.txt");
    ofile << *computeVelocity(lattice);
}


int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    //defaultMultiBlockPolicy3D().toggleBlockingCommunication(true);

    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 10.,   // Re
            100,        // N
            1.,        // lx
            1.,        // ly
            1.         // lz
    );
    const T logT     = (T)1/(T)10;
    const T imSave   = (T)1/(T)10;
    const T vtkSave  = (T)1/(T)100;
    const T maxT     = (T)10.1;

    pcout << "omega= " << parameters.getOmega() << endl;
    writeLogFile(parameters, "3D diagonal cavity");

    T omega = parameters.getOmega();

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(), parameters.getNz(),
            new BGKdynamics<T,DESCRIPTOR>(omega) );


    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        //= createInterpBoundaryCondition3D<T,DESCRIPTOR>();
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    cavitySetup(lattice, parameters, *boundaryCondition);

    // Transient statistics are statistics like "mean value" and "standard deviation" that
    // refer to time. They must be registered, initialized and updated properly during the
    // course of the simulation.
    // Here we define an object of type TransientStatistics3D as a manager for all transient
    // statistics operations. This object is bound with a specific lattice and domain. If one
    // wants statistics over more domains, she should instantiate more objects appropriately.
    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();
    Box3D domainForTransientStatistics(nx/2-1, nx/2+1, 0, ny-1, 0, nz-1);
    TransientStatistics3D<T,DESCRIPTOR> transientStatistics(lattice, domainForTransientStatistics);
    initializeTransientStatistics(transientStatistics);

    T previousIterationTime = T();
    // Loop over main time iteration.
    for (plint iT=0; iT<parameters.nStep(maxT); ++iT) {
        global::timer("mainLoop").restart();

        if (iT%parameters.nStep(imSave)==0) {
            pcout << "Writing Gif ..." << endl;
            writeGifs(lattice, parameters, iT);
        }

        if (iT%parameters.nStep(vtkSave)==0 && iT>0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(lattice, parameters, iT);

            // Next we update and save all transient statistics.
            transientStatistics.update();
            string path("");
            string domainName("slice_x");
            plint namePadding = 6;
            Array<T,3> physicalLocation; physicalLocation.resetToZero();
            T rhoFluid = 1.0;
            T pressureOffset = 1.0e5;
            T rhoLB = 1.0;
            transientStatistics.output(path, domainName, iT, namePadding, parameters.getDeltaX(), parameters.getDeltaT(),
                    physicalLocation, rhoFluid, pressureOffset, rhoLB);
        }

        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; t=" << iT*parameters.getDeltaT();
        }

        // Execute a time iteration.
        lattice.collideAndStream();

        // Access averages from internal statistics ( their value is defined
        //   only after the call to lattice.collideAndStream() )
        if (iT%parameters.nStep(logT)==0) {
            pcout << "; av energy="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho="
                  << setprecision(10) << getStoredAverageDensity<T>(lattice) << endl;
            pcout << "Time spent during previous iteration: "
                  << previousIterationTime << endl;
        }

        previousIterationTime = global::timer("mainLoop").stop();
    }

    delete boundaryCondition;
}

