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

/* Code 2.2 in the Palabos tutorial
 */

#include "palabos2D.h"
#include "palabos2D.hh"
#include <vector>
#include <iostream>
#include <iomanip>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor

/// Describe the geometry of the half-circular channel, used in tutorial 2.
template<typename T>
class BounceBackNodes : public DomainFunctional2D {
public:
    BounceBackNodes(plint N, plint radius)
        : cx(N/2),
          cy(N/2),
          innerR(radius),
          outerR(N/2)
    { }
    /// Return true for all cells outside the channel, on which bounce-back
    ///  dynamics must be instantiated.
    virtual bool operator() (plint iX, plint iY) const {
        T rSqr = util::sqr(iX-cx) + util::sqr(iY-cy);
        return rSqr <= innerR*innerR || rSqr >= outerR*outerR;
    }
    virtual BounceBackNodes<T>* clone() const {
        return new BounceBackNodes<T>(*this);
    }
private:
    plint cx;      //< X-position of the center of the half-circle.
    plint cy;      //< Y-position of the center of the half-circle.
    plint innerR;  //< Outer radius of the half-circle.
    plint outerR;  //< Inner radius of the half-circle.
};


void halfCircleSetup (
        MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint N, plint radius,
        OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )
{
    // The channel is pressure-driven, with a difference deltaRho
    //   between inlet and outlet.
    T deltaRho = 1.e-2;
    T rhoIn  = 1. + deltaRho/2.;
    T rhoOut = 1. - deltaRho/2.;

    Box2D inlet (0,     N/2, N/2, N/2);
    Box2D outlet(N/2+1, N,   N/2, N/2);

    boundaryCondition.addPressureBoundary1P(inlet, lattice);
    boundaryCondition.addPressureBoundary1P(outlet, lattice);

    // Specify the inlet and outlet density.
    setBoundaryDensity (lattice, inlet, rhoIn);
    setBoundaryDensity (lattice, outlet, rhoOut);

    // Create the initial condition.
    Array<T,2> zeroVelocity((T)0.,(T)0.);
    T constantDensity = (T)1;
    initializeAtEquilibrium (
       lattice, lattice.getBoundingBox(), constantDensity, zeroVelocity );

    defineDynamics(lattice, lattice.getBoundingBox(),
                   new BounceBackNodes<T>(N, radius),
                   new BounceBack<T,DESCRIPTOR>);

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

    // Parameters of the simulation
    plint N         = 400;    // Use a 400x200 domain.
    plint maxT      = 20001;
    plint imageIter = 1000;
    T omega        = 1.;
    plint radius    = N/3;    // Inner radius of the half-circle.

    // Parameters for the creation of the multi-block.
    
    // d is the width of the block which is exempted from the full domain.
    plint d = (plint) (2.*std::sqrt((T)util::sqr(radius)-(T)util::sqr(N/4.)));
    plint x0 = (N-d)/2 + 1;  // Begin of the exempted block.
    plint x1 = (N+d)/2 - 1;  // End of the exempted block.

    // Create a block distribution with the three added blocks.
    plint envelopeWidth = 1;
    SparseBlockStructure2D sparseBlock(N+1, N/2+1);
    sparseBlock.addBlock(Box2D(0, x0,      0, N/2),   sparseBlock.nextIncrementalId());
    sparseBlock.addBlock(Box2D(x0+1, x1-1, 0, N/4+1), sparseBlock.nextIncrementalId());
    sparseBlock.addBlock(Box2D(x1, N,      0, N/2),   sparseBlock.nextIncrementalId());

    // Instantiate the multi-block, based on the created block distribution and
    // on default parameters.
    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
        MultiBlockManagement2D (
            sparseBlock,
            defaultMultiBlockPolicy2D().getThreadAttribution(), envelopeWidth ),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiCellAccess<T,DESCRIPTOR>(),
        new BGKdynamics<T,DESCRIPTOR>(omega)
    );

    pcout << getMultiBlockInfo(lattice) << std::endl;

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    halfCircleSetup(lattice, N, radius, *boundaryCondition);

    // Main loop over time iterations.
    for (plint iT=0; iT<maxT; ++iT) {
        if (iT%imageIter==0) {
            pcout << "Saving Gif at time step " << iT << endl;
            writeGifs(lattice, iT);
        }
        lattice.collideAndStream();
    }

    delete boundaryCondition;
}
