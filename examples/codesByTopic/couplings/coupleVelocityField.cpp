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
#include "poiseuille.h"
#include "poiseuille.hh"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D2Q9Descriptor

template<typename T>
class SetToPoiseuilleVelocityFunctional : public BoxProcessingFunctional2D_T<T,2> {
public:
    SetToPoiseuilleVelocityFunctional(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    virtual void process(Box2D domain, TensorField2D<T,2>& velocity)
    {
        Dot2D offset = velocity.getLocation();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint absoluteY = iY + offset.y;
                T ux = poiseuilleVelocity(absoluteY,parameters);
                velocity.get(iX,iY) = Array<T,2>(ux, T());
            }
        }
    }
    virtual SetToPoiseuilleVelocityFunctional<T>* clone() const {
        return new SetToPoiseuilleVelocityFunctional(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
private:
    IncomprFlowParam<T> parameters;
};

void setToPoiseuilleProfile(MultiTensorField2D<T,2>& velocity,
                            IncomprFlowParam<T> parameters)
{
    applyProcessingFunctional(new SetToPoiseuilleVelocityFunctional<T>(parameters),
                              velocity.getBoundingBox(), velocity);
}

template<typename T, template<typename U> class Descriptor>
class BoundaryFromVelocityFunctional2D : public BoxProcessingFunctional2D_LT<T,Descriptor,T,2> {
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                         TensorField2D<T,2>& velocity)
    {
        Dot2D displacement = computeRelativeDisplacement(lattice, velocity);
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                lattice.get(iX,iY).defineVelocity(velocity.get(iX+displacement.x,iY+displacement.y));
            }
        }
    }
    virtual BoundaryFromVelocityFunctional2D<T,Descriptor>* clone() const {
        return new BoundaryFromVelocityFunctional2D(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::dataStructure;
        modified[1] = modif::nothing;
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
};

void createBoundariesFromVelocityField(MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                                       MultiTensorField2D<T,2>& velocity)
{
    applyProcessingFunctional( new BoundaryFromVelocityFunctional2D<T,DESCRIPTOR>,
                               lattice.getBoundingBox(),
                               lattice, velocity );
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters (
        (T) 1e-2,  // uMax
        (T) 10.,   // Re
        30,        // N
        2.,        // lx
        1.         // ly 
    );

    plint nx = parameters.getNx();
    plint ny = parameters.getNy();

    writeLogFile(parameters, "Poiseuille flow");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
              nx, ny,
              new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    MultiTensorField2D<T,2> velocity(nx,ny);

    setToPoiseuilleProfile(velocity, parameters);

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();
    boundaryCondition -> setVelocityConditionOnBlockBoundaries(lattice);
    pcout << "start" << endl;
    createBoundariesFromVelocityField(lattice, velocity);
    pcout << "stop" << endl;
    lattice.initialize();

    // The following command opens a text-file, in which the velocity-profile
    // in the middle of the channel will be written and several successive
    // time steps. Note the use of plb_ofstream instead of the standard C++
    // ofstream, which is required to guarantee a consistent behavior in MPI-
    // parallel programs.
    plb_ofstream successiveProfiles("velocityProfiles.dat");

    // Main loop over time steps.
    for (plint iT=0; iT<10000; ++iT) {
        if (iT%1000==0) {

            ImageWriter<T> imageWriter("leeloo");
            imageWriter.writeScaledGif(createFileName("uSqr", iT, 6),
                                       *computeKineticEnergy(lattice),
                                       600, 600 );
            pcout << "At iteration step " << iT
                  << ", the density along the channel is " << endl;
            pcout << setprecision(7)
                  << *computeDensity(lattice, Box2D(0, nx-1, ny/2, ny/2))
                  << endl << endl;

            Box2D profileSection(nx/2, nx/2, 0, ny-1);
            successiveProfiles
                << setprecision(4)
                  // (2) Convert from lattice to physical units.
                << *multiply (
                       parameters.getDeltaX() / parameters.getDeltaT(),
                  // (1) Compute velocity norm along the chosen section.
                       *computeVelocityNorm (lattice, profileSection) )
                << endl;

        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();
    }

    delete boundaryCondition;
}
