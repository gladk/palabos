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
#include "complexDynamics/smagorinskyDynamics3D.h"

#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "io/parallelIO.h"

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

void randomIniCondition(plint iX, plint iY, plint iZ, T randomValue, T& rho, Array<T,3>& velocity) {
    velocity.resetToZero();
    rho = (T)1 + 1.e-2*randomValue;
}

void cavitySetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                  IncomprFlowParam<T> const& parameters,
                  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    Box3D topLid = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);
    Box3D everythingButTopLid = Box3D(0, nx-1, 0, ny-2, 0, nz-1);

    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    T u = std::sqrt((T)2)/(T)2 * parameters.getLatticeU();
    initializeAtRandomEquilibrium(lattice, everythingButTopLid, randomIniCondition);
    initializeAtEquilibrium(lattice, topLid, (T)1., Array<T,3>(u,(T)0.,u) );
    setBoundaryVelocity(lattice, topLid, Array<T,3>(u,(T)0.,u) );

    lattice.initialize();
}

void writeGifs(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
               IncomprFlowParam<T> const& parameters, int iter)
{
    const plint imSize = 600;
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();

    Box3D zSlice(0, nx-1, 0, ny-1, nz/2, nz/2);
    Box3D ySlice(0, nx-1, ny/2, ny/2, 0, nz-1);
    Box3D extendedZslice(0, nx-1, 0, ny-1, nz/2-1, nz/2+1);
    ImageWriter<T> imageWriter("leeloo");

    imageWriter.writeScaledGif( createFileName("zSlice_u", iter, 6),
                                *computeVelocityNorm (lattice, zSlice),
                                imSize, imSize );
    imageWriter.writeScaledGif( createFileName("ySlice_u", iter, 6),
                                *computeVelocityNorm (lattice, ySlice),
                                imSize, imSize );
    imageWriter.writeScaledGif( createFileName("zSlice_omega", iter, 6),
                                *computeNorm (
                                    *computeVorticity(*computeVelocity(lattice, extendedZslice)),
                                    zSlice ),
                                imSize, imSize );

}

void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
}


int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 1.e4,  // Re
            90,        // N
            1.,        // lx
            1.,        // ly
            1.         // lz
    );
    const T logT     = (T)1/(T)500;
    const T imSave   = (T)1/(T)10;
    const T vtkSave  = (T)10;
    const T maxT     = (T)50;

    writeLogFile(parameters, "3D diagonal cavity");

    T cSmago = 0.14;
    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(), parameters.getNz(),
            //new SmagorinskyBGKdynamics<T,DESCRIPTOR>(parameters.getOmega(), cSmago) );
            new SmagorinskyRegularizedDynamics<T,DESCRIPTOR>(parameters.getOmega(), cSmago) );
            //new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    // Uncomment the following line to instantiate the Smagorinsky LES model,
    //   if the background-dynamics (i.e. the dynamics given to "lattice" in
    //   the constructor) is BGKdynamics instead of SmagorinskyBGKdynamics.
    
      // instantiateStaticSmagorinsky(lattice, lattice.getBoundingBox(), cSmago);

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    cavitySetup(lattice, parameters, *boundaryCondition);

    T previousIterationTime = T();

    // Main loop over time iterations.
    for (plint iT=0; iT<parameters.nStep(maxT); ++iT) {
        global::timer("mainLoop").restart();

        if (iT%parameters.nStep(imSave)==0) {
            pcout << "Writing Gif ..." << endl;
            writeGifs(lattice, parameters, iT);
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
