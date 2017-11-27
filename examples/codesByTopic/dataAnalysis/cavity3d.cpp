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
#include "palabos3D.hh"   // Use full generic version

#include <vector>
#include <cmath>
#include <cstdlib>
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

    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    T u = std::sqrt((T)2)/(T)2 * parameters.getLatticeU();
    initializeAtEquilibrium(lattice, everythingButTopLid, (T) 1., Array<T,3>((T)0.,(T)0.,(T)0.) );
    initializeAtEquilibrium(lattice, topLid, (T) 1., Array<T,3>(u,(T)0.,u) );
    setBoundaryVelocity(lattice, topLid, Array<T,3>(u,(T)0.,u) );

    lattice.initialize();
}

void analyzeStrainRate(MultiBlockLattice3D<T,DESCRIPTOR>& lattice)
{
    Box3D slice(0, lattice.getNx()-1, 0, lattice.getNy()-1, lattice.getNz()/2, lattice.getNz()/2);
    auto_ptr<MultiTensorField3D<T,6> > strainRate_LB (
            computeStrainRateFromStress(lattice) );
    auto_ptr<MultiTensorField3D<T,6> > strainRate_FD (
            computeStrainRate(*computeVelocity(lattice)) );
    auto_ptr<MultiScalarField3D<T> > Snorm_LB (
            computeSymmetricTensorNorm(*strainRate_LB) );
    auto_ptr<MultiScalarField3D<T> > Snorm_FD (
            computeSymmetricTensorNorm(*strainRate_FD) );

    const plint imSize = 400;
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif( "S_LB",
            *extractSubDomain(*Snorm_LB, slice), imSize,imSize );
    imageWriter.writeScaledGif( "S_FD",
            *extractSubDomain(*Snorm_FD, slice), imSize,imSize );
    imageWriter.writeScaledGif( "S_LB_Trace",
            *extractSubDomain(*computeSymmetricTensorTrace(*strainRate_LB), slice), imSize,imSize );
    if (global::mpi().isMainProcessor()) {
        int err = 0;
        err = system("convert +append tmp/S_FD.gif tmp/S_LB.gif tmp/S_LB_Trace.gif tmp/strain.gif");
        if (err != 0)
            exit(err);
        err = system("/bin/rm tmp/S_FD.gif tmp/S_LB.gif tmp/S_LB_Trace.gif");
        if (err != 0)
            exit(err);
    }

    auto_ptr<MultiScalarField3D<T> > differenceSqr (
            computeSymmetricTensorNorm (
                *subtract (
                    *strainRate_LB, *strainRate_FD ) ) );

    imageWriter.writeScaledGif( "differenceSqr",
            *extractSubDomain(*differenceSqr, slice), imSize,imSize );

    T maxSnorm_LB = computeMax(*Snorm_LB);
    T maxSnorm_FD = computeMax(*Snorm_FD);
    pcout << "Max LB " << maxSnorm_LB << endl;
    pcout << "Max FD " << maxSnorm_FD << endl;
    T differenceMax = std::sqrt(computeMax(*differenceSqr));
    T differenceMin = std::sqrt(computeMin(*differenceSqr));

    pcout << "Absolute difference range: ["
          << differenceMin << "," << differenceMax << "]" << endl;
    pcout << "Relative difference range: ["
          << differenceMin/maxSnorm_LB << ","
          << differenceMax/maxSnorm_LB << "]" << endl;
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters(
            (T) 2e-2,  // uMax
            (T) 40.,   // Re
            40,        // N
            1.,        // lx
            1.,        // ly
            1.         // lz
    );
    const T logT     = (T)1/(T)100;
    const T imSave   = (T)0.2;
    const T maxT     = (T)10.1;

    writeLogFile(parameters, "3D diagonal cavity");

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(), parameters.getNz(),
            new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        //= createInterpBoundaryCondition3D<T,DESCRIPTOR>();
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    cavitySetup(lattice, parameters, *boundaryCondition);

    // Main loop over time iterations.
    for (plint iT=0; iT<parameters.nStep(maxT); ++iT) {
        if (iT%parameters.nStep(imSave)==0) {
            pcout << "Processing analysis of the strain rate ..." << endl
                  << "=========================================="
                  << endl << endl;
            analyzeStrainRate(lattice);
        }

        if (iT%parameters.nStep(logT)==0) {
            pcout << "Compute averages at time step " << iT << endl;
            pcout << "================================" << endl;
            pcout << "t=" << iT*parameters.getDeltaT() << endl;
            // Compute averages from the lattice. This takes a certain amount
            //   of time to be computed, because all lattice cells must be
            //   accessed.
            pcout << "Av energy="
                  << setprecision(10) << computeAverageEnergy(lattice)
                  << "; av rho="
                  << setprecision(10) << computeAverageDensity(lattice)
                  << endl << endl;
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();

        if (iT%parameters.nStep(logT)==0) {
            pcout << "Access internal statistics for time step " << iT << endl;
            pcout << "===========================================" << endl;
            // Display averages by accessing internal statistics of the lattice.
            //   This takes practically no time, because the averages where
            //   computed automatically during the lattice iteration cycle.
            pcout << "Av energy="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho="
                  << setprecision(10) << getStoredAverageDensity<T>(lattice)
                  << endl << endl;
        }

    }

    delete boundaryCondition;
}
