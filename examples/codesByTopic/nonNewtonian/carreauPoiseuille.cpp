/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2012 FlowKit Sarl
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
#include "poiseuille.h"
#include "poiseuille.hh"
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

void initialAndBoundaryCondition( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                                  IncomprFlowParam<T> const& parameters,
                                  OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();

    // Velocity boundary condition on bottom wall. 
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, nx-1, 0, 0) );
    // Velocity boundary condition on top wall. 
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, nx-1, ny-1, ny-1) );

    // Pressure condition on inlet.
    boundaryCondition.setPressureConditionOnBlockBoundaries (
            lattice, Box2D(0,0, 1,ny-2) );
    // Pressure condition on outlet.
    boundaryCondition.setPressureConditionOnBlockBoundaries (
            lattice, Box2D(nx-1,nx-1, 1,ny-2) );

    setBoundaryDensity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleDensity<T,DESCRIPTOR>(parameters) );
    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<T>(parameters) );
    initializeAtEquilibrium (
           lattice, lattice.getBoundingBox(),
           PoiseuilleVelocityAndDensity<T,DESCRIPTOR>(parameters) );

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
            5.,        // lx
            1.         // ly 
    );
    const T logT     = (T)0.02;
    const T imSave   = (T)0.06;
    const T maxT     = (T)3.1;

    writeLogFile(parameters, "Poiseuille flow");

    global::CarreauParameters().setNu0(parameters.getLatticeNu());
    global::CarreauParameters().setLambda(1.);
    global::CarreauParameters().setExponent(0.3);

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
              parameters.getNx(), parameters.getNy(),
              new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        //boundaryCondition = createInterpBoundaryCondition2D<T,DESCRIPTOR>();
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    initialAndBoundaryCondition(lattice, parameters, *boundaryCondition);

    setCompositeDynamics (
            lattice,
            lattice.getBoundingBox(),
            new CarreauDynamics<T,DESCRIPTOR,0>(new NoDynamics<T,DESCRIPTOR>) );
            
    // Main loop over time iterations.
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {
        if (iT%parameters.nStep(imSave)==0 && iT>0) {
            pcout << "Saving Gif ..." << endl;
            writeGifs(lattice, iT);
        }

        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; lattice time=" << lattice.getTimeCounter().getTime()
                  << "; t=" << iT*parameters.getDeltaT();
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();

        if (iT%parameters.nStep(logT)==0) {
            pcout << "; av energy="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho="
                  << getStoredAverageDensity<T>(lattice) << endl;
        }
    }

    delete boundaryCondition;
}
