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
#include "palabos2D.h"
#include "palabos2D.hh"
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace plb::descriptors;
using namespace std;

plint nx = 40;
plint ny = 40;
plint nz = 40;
double omega = 1.;

Box3D domain1_3d(5,15,5,15,5,15);
Box3D domain2_3d(20,30,20,30,20,30);

Box2D domain1_2d(5,15,5,15);
Box2D domain2_2d(20,30,20,30);


void copyBlockLattice3D() {
    MultiBlockLattice3D<double, D3Q19Descriptor> lattice3d_1 (
            nx,ny,nz, new BGKdynamics<double,D3Q19Descriptor>(omega) );
    MultiBlockLattice3D<double, D3Q19Descriptor> lattice3d_2 (
            nx,ny,nz, new BGKdynamics<double,D3Q19Descriptor>(omega) );

    lattice3d_1.collideAndStream();

    // 1. SAME-POSITION COPIES.
    // Copy all populations and external scalars.
    copyPopulations(lattice3d_1, lattice3d_2, domain1_3d);
    // Copy all populations, external scalars, and content of dynamics objects.
    copyAll(lattice3d_1, lattice3d_2, domain1_3d);
    defineDynamics(lattice3d_1, domain1_3d,
                   new BounceBack<double,D3Q19Descriptor>(1.));
    // Copy everything, and recreate a dynamics object at the destination.
    copyRegenerate(lattice3d_1, lattice3d_2, domain1_3d);

    // 2. CROSS-POSITION COPIES.
    // Copy all populations and external scalars.
    copyPopulations(lattice3d_1, domain1_3d, lattice3d_2, domain2_3d);
    // Copy all populations, external scalars, and content of dynamics objects.
    copyAll(lattice3d_1, domain1_3d, lattice3d_2, domain2_3d);
    // Copy everything, and recreate a dynamics object at the destination.
    copyRegenerate(lattice3d_1, domain1_3d, lattice3d_2, domain2_3d);

    // 3. CREATING CLONES.
    // Create a reduced-size clone, correponding to a sub-domain of the original lattice.
    std::auto_ptr<MultiBlockLattice3D<double,D3Q19Descriptor> >
        lattice3d_3 = clone(lattice3d_1,domain1_3d);
}

void copyBlockLattice2D() {
    MultiBlockLattice2D<double, D2Q9Descriptor> lattice2d_1 (
            nx,ny, new BGKdynamics<double,D2Q9Descriptor>(omega) );
    MultiBlockLattice2D<double, D2Q9Descriptor> lattice2d_2 (
            nx,ny, new BGKdynamics<double,D2Q9Descriptor>(omega) );

    lattice2d_1.collideAndStream();

    // 1. SAME-POSITION COPIES.
    // Copy all populations and external scalars.
    copyPopulations(lattice2d_1, lattice2d_2, domain1_2d);
    // Copy all populations, external scalars, and content of dynamics objects.
    copyAll(lattice2d_1, lattice2d_2, domain1_2d);
    defineDynamics(lattice2d_1, domain1_2d,
                   new BounceBack<double,D2Q9Descriptor>(1.));
    // Copy everything, and recreate a dynamics object at the destination.
    copyRegenerate(lattice2d_1, lattice2d_2, domain1_2d);

    // 2. CROSS-POSITION COPIES.
    // Copy all populations and external scalars.
    copyPopulations(lattice2d_1, domain1_2d, lattice2d_2, domain2_2d);
    // Copy all populations, external scalars, and content of dynamics objects.
    copyAll(lattice2d_1, domain1_2d, lattice2d_2, domain2_2d);
    // Copy everything, and recreate a dynamics object at the destination.
    copyRegenerate(lattice2d_1, domain1_2d, lattice2d_2, domain2_2d);

    // 3. CREATING CLONES.
    // Create a reduced-size clone, correponding to a sub-domain of the original lattice.
    std::auto_ptr<MultiBlockLattice2D<double,D2Q9Descriptor> >
        lattice2d_3 = clone(lattice2d_1,domain1_2d);
}

void copyScalarField3D() {
    MultiScalarField3D<double> scalar3d_1(nx,ny,nz);
    MultiScalarField3D<double> scalar3d_2(nx,ny,nz);
    MultiScalarField3D<int> scalar3d_3(nx,ny,nz);
    setToConstant(scalar3d_1, scalar3d_1.getBoundingBox(), 1.);
    setToConstant(scalar3d_2, scalar3d_2.getBoundingBox(), 2.);
    setToConstant(scalar3d_3, scalar3d_3.getBoundingBox(), 3);

    // 1. SAME-POSITION COPIES.
    plb::copy(scalar3d_1, scalar3d_2, scalar3d_1.getBoundingBox() );
    // Copy with type conversion.
    plb::copy(scalar3d_2, scalar3d_3, scalar3d_2.getBoundingBox() );
    // 2. CROSS-POSITION COPIES.
    plb::copy(scalar3d_1, domain1_3d, scalar3d_2, domain2_3d);
    // 3. CREATING CLONES.
    std::auto_ptr<MultiScalarField3D<double> >
        scalar3d_4 = clone(scalar3d_2,domain1_3d);
}

void copyScalarField2D() {
    MultiScalarField2D<double> scalar2d_1(nx,ny);
    MultiScalarField2D<double> scalar2d_2(nx,ny);
    MultiScalarField2D<int> scalar2d_3(nx,ny);
    setToConstant(scalar2d_1, scalar2d_1.getBoundingBox(), 1.);
    setToConstant(scalar2d_2, scalar2d_2.getBoundingBox(), 2.);
    setToConstant(scalar2d_3, scalar2d_3.getBoundingBox(), 2);

    // 1. SAME-POSITION COPIES.
    plb::copy(scalar2d_1, scalar2d_2, scalar2d_1.getBoundingBox() );
    // Copy with type conversion.
    plb::copy(scalar2d_2, scalar2d_3, scalar2d_2.getBoundingBox() );
    // 2. CROSS-POSITION COPIES.
    plb::copy(scalar2d_1, domain1_2d, scalar2d_2, domain2_2d);
    // 3. CREATING CLONES.
    std::auto_ptr<MultiScalarField2D<int> >
        scalar2d_4 = clone(scalar2d_3,domain1_2d);
}


int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    copyBlockLattice3D();
    copyScalarField3D();
}
