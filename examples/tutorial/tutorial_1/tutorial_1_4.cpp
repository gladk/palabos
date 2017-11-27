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

/* Code 1.4 in the Palabos tutorial
 */

#include "palabos2D.h"
#include "palabos2D.hh"
#include <iostream>
#include <iomanip>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor

const plint maxIter = 100;  // Iterate during 100 steps.
const plint nx = 600;       // Choice of lattice dimensions.
const plint ny = 600;
const T omega = 1.;        // Choice of the relaxation parameter

T rho0 = 1.; // All cells have initially density rho ...
// .. except for those in the box "centralSquare" which have density
//    rho+deltaRho
T deltaRho = 1.e-4;
Array<T,2> u0(0,0);

void initializeRhoOnCircle(plint iX, plint iY, T& rho, Array<T,2>& u) {
    plint radius = nx/6;
    plint centerX = nx/3;
    plint centerY = ny/4;
    u = u0;
    if( (iX-centerX)*(iX-centerX) + (iY-centerY)*(iY-centerY) < radius*radius) {
        rho = rho0 + deltaRho;
    }
    else {
        rho = rho0;
    }
}

// Initialize the lattice at zero velocity and constant density, except
//   for a slight density excess on a square sub-domain.
void defineInitialDensityAtCenter(MultiBlockLattice2D<T,DESCRIPTOR>& lattice)
{
    // Initialize constant density everywhere.
    initializeAtEquilibrium (
           lattice, lattice.getBoundingBox(), rho0, u0 );

    // And slightly higher density in the central box.
    initializeAtEquilibrium (
           lattice, lattice.getBoundingBox(), initializeRhoOnCircle );

    lattice.initialize();
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
           nx, ny, new BGKdynamics<T,DESCRIPTOR>(omega) );

    lattice.periodicity().toggleAll(true); // Set periodic boundaries.

    defineInitialDensityAtCenter(lattice);

    // First part: loop over time iterations.
    for (plint iT=0; iT<maxIter; ++iT) {
        lattice.collideAndStream();
    }

    // Second part: Data analysis.
    Array<T,2> velocity;
    lattice.get(nx/2, ny/2).computeVelocity(velocity);
    pcout << "Velocity in the middle of the lattice: ("
          << velocity[0] << "," << velocity[1] << ")" << endl;

    pcout << "Velocity norm along a horizontal line: " << endl;
    Box2D line(0, 100, ny/2, ny/2);
    pcout << setprecision(3) << *computeVelocityNorm(*extractSubDomain(lattice, line)) << endl;

    plb_ofstream ofile("profile.dat");
    ofile << setprecision(3) << *computeVelocityNorm(*extractSubDomain(lattice, line)) << endl;

    pcout << "Average density in the domain: " << computeAverageDensity(lattice) << endl;
    pcout << "Average energy in the domain: " << computeAverageEnergy(lattice) << endl;
}
