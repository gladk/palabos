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
#include "poiseuille.h"
#include "poiseuille.hh"
#include "cylinder.h"
#include "cylinder.hh"

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;

// Uncomment the two following lines for BGK dynamics
    //#define DESCRIPTOR D2Q9Descriptor
    //typedef BGKdynamics<T,DESCRIPTOR> BackgroundDynamics;

// Uncomment the two following lines for so-called incompressible BGK dynamics
    //#define DESCRIPTOR D2Q9Descriptor
    //typedef IncBGKdynamics<T,DESCRIPTOR> BackgroundDynamics;

// Uncomment the two following lines for Regularized BGK dynamics
    //#define DESCRIPTOR D2Q9Descriptor
    //typedef RegularizedBGKdynamics<T,DESCRIPTOR> BackgroundDynamics;

// Uncomment the two following lines for MRT dynamics
     #define DESCRIPTOR MRTD2Q9Descriptor
     typedef MRTdynamics<T,DESCRIPTOR> BackgroundDynamics;

// Uncomment the two following lines for Entropic dynamics
    //#define DESCRIPTOR D2Q9Descriptor
    //typedef EntropicDynamics<T,DESCRIPTOR> BackgroundDynamics;


void defineCylinderGeometry( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                             IncomprFlowParam<T> const& parameters )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();

    plint cx      = nx/4;
    plint cy      = ny/2+ny/10;
    plint radius  = cy/4;

    createCylinder(lattice, cx, cy, radius);
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
            (T) 300.,  // Re
            100,       // N
            5.,        // lx
            1.         // ly 
    );
    const T logT     = (T)0.02;
    const T imSave   = (T)0.1;
    const T maxT     = (T)10.1;

    writeLogFile(parameters, "Poiseuille flow");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(),
            new BackgroundDynamics(parameters.getOmega()) );

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    defineCylinderGeometry(lattice, parameters);

    createPoiseuilleBoundaries(lattice, parameters, *boundaryCondition);
    createPoiseuilleInitialValues(lattice, parameters);

    // Main loop over time iterations.
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {
        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; lattice time=" << lattice.getTimeCounter().getTime()
                  << "; t=" << iT*parameters.getDeltaT()
                  << "; av energy="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho="
                  << getStoredAverageDensity<T>(lattice) << endl;
        }

        if (iT%parameters.nStep(imSave)==0) {
            pcout << "Saving Gif ..." << endl;
            writeGifs(lattice, iT);
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();
    }

    delete boundaryCondition;
}
