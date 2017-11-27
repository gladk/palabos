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
 * Simulation of a 3D Rayleigh-Taylor instability, which describes a symmetry breakdown,
 * as a heavy fluid, initially located on top of a light one, starts penetrating the
 * latter. This is an illutration of the Shan/Chen multi-component model. Note that the
 * single-component multi-phase model is also implemented in Palabos, and a corresponding
 * sample program is provided. Also, note that the multi-component model can be used to
 * couple more than two phases, using the same approach as the one shown here.
 **/

#include "palabos3D.h"
#include "palabos3D.hh"
#include <cstdlib>
#include <iostream>

using namespace plb;
using namespace std;

// Use double-precision arithmetics
typedef double T;
// Use a grid which additionally to the f's stores three variables for
//   the external force term.
#define DESCRIPTOR descriptors::ForcedShanChenD3Q19Descriptor

/// Initial condition: heavy fluid on top, light fluid on bottom.
/** This functional is going to be used as an argument to the function "applyIndexed",
 *  to setup the initial condition. For efficiency reasons, this approach should
 *  always be preferred over explicit space loops in end-user codes.
 */
template<typename T, template<typename U> class Descriptor>
class TwoLayerInitializer : public OneCellIndexedWithRandFunctional3D<T,Descriptor> {
public:
    TwoLayerInitializer(plint ny_, bool topLayer_)
        : ny(ny_),
          topLayer(topLayer_)
    { }
    TwoLayerInitializer<T,Descriptor>* clone() const {
        return new TwoLayerInitializer<T,Descriptor>(*this);
    }
    virtual void execute(plint iX, plint iY, plint iZ, T rand_val, Cell<T,Descriptor>& cell) const {
        T densityFluctuations = 1.e-2;
        T almostNoFluid       = 1.e-4;
        Array<T,3> zeroVelocity (0.,0.,0.);

        T rho = (T)1;
        // Add a random perturbation to the initial condition to instantiate the
        //   instability.
        if ( (topLayer && iY>ny/2) || (!topLayer && iY <= ny/2) ) {
            rho += rand_val * densityFluctuations;
        }
        else {
            rho = almostNoFluid;
        }

        iniCellAtEquilibrium(cell, rho, zeroVelocity);
    }
private:
    plint ny;
    bool topLayer;
};

void rayleighTaylorSetup( MultiBlockLattice3D<T, DESCRIPTOR>& heavyFluid,
                          MultiBlockLattice3D<T, DESCRIPTOR>& lightFluid,
                          T rho0, T rho1,
                          T force )
{
    // The setup is: periodicity along horizontal direction, bounce-back on top
    // and bottom. The upper half is initially filled with fluid 1 + random noise,
    // and the lower half with fluid 2. Only fluid 1 experiences a forces,
    // directed downwards.
    plint nx = heavyFluid.getNx();
    plint ny = heavyFluid.getNy();
    plint nz = heavyFluid.getNz();
    
    // Bounce-back on bottom wall (where the light fluid is, initially).
    defineDynamics(heavyFluid, Box3D(0,nx-1, 0,0, 0,nz-1), new BounceBack<T, DESCRIPTOR>(rho0) );
    defineDynamics(lightFluid, Box3D(0,nx-1, 0,0, 0,nz-1), new BounceBack<T, DESCRIPTOR>(rho1) );
    // Bounce-back on top wall (where the heavy fluid is, initially).
    defineDynamics(heavyFluid, Box3D(0,nx-1, ny-1,ny-1, 0,nz-1), new BounceBack<T, DESCRIPTOR>(rho1) );
    defineDynamics(lightFluid, Box3D(0,nx-1, ny-1,ny-1, 0,nz-1), new BounceBack<T, DESCRIPTOR>(rho0) );
   
    // Initialize top layer.
    applyIndexed(heavyFluid, Box3D(0, nx-1, 0, ny-1, 0, nz-1),
                 new TwoLayerInitializer<T,DESCRIPTOR>(ny, true) );
    // Initialize bottom layer.
    applyIndexed(lightFluid, Box3D(0, nx-1, 0, ny-1, 0, nz-1),
                 new TwoLayerInitializer<T,DESCRIPTOR>(ny, false) );

    // Let's have gravity acting on the heavy fluid only. This represents a situation
    //   where the molecular mass of the light species is very small, and thus the
    //   action of gravity on this species is negligible.
    setExternalVector(heavyFluid, heavyFluid.getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T,3>(0.,-force,0.));
    setExternalVector(lightFluid, lightFluid.getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T,3>(0.,0.,0.));

    lightFluid.initialize();
    heavyFluid.initialize();
}

void writePpms(MultiBlockLattice3D<T, DESCRIPTOR>& heavyFluid,
               MultiBlockLattice3D<T, DESCRIPTOR>& lightFluid, plint iT)
{
    const plint nx = heavyFluid.getNx();
    const plint ny = heavyFluid.getNy();
    const plint nz = heavyFluid.getNz();
    Box3D slice(0, nx-1, 0, ny-1, nz/2, nz/2);

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledPpm(createFileName("rho_heavy_", iT, 6),
                               *computeDensity(heavyFluid, slice));
    imageWriter.writeScaledPpm(createFileName("rho_light_", iT, 6),
                               *computeDensity(lightFluid, slice));
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    srand(global::mpi().getRank());
    
    const T omega1 = 1.0;
    const T omega2 = 1.0;
    const plint nx   = 225;
    const plint ny   = 75;
    const plint nz   = 225;
    const T G      = 1.2;
    T force        = 0.15/(T)ny; //1.e-3;
    const plint maxIter  = 10000;
    const plint saveIter = 100;
    const plint statIter = 50;

    // Use regularized BGK dynamics to improve numerical stability (but note that
    //   BGK dynamics works well too).
    MultiBlockLattice3D<T, DESCRIPTOR> heavyFluid (
            nx,ny,nz, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega1) );
    MultiBlockLattice3D<T, DESCRIPTOR> lightFluid (
            nx,ny,nz, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega2) );

    // Make x- and z-directions periodic.
    heavyFluid.periodicity().toggle(0, true);
    heavyFluid.periodicity().toggle(2, true);
    lightFluid.periodicity().toggle(0, true);
    lightFluid.periodicity().toggle(2, true);

    T rho1 = 1.; // Fictitious density experienced by the partner fluid on a Bounce-Back node.
    T rho0 = 0.; // Fictitious density experienced by the partner fluid on a Bounce-Back node.
    
    // Store a pointer to all lattices (two in the present application) in a vector to
    //   create the Shan/Chen coupling therm. The heavy fluid being at the first place
    //   in the vector, the coupling term is going to be executed at the end of the call
    //   to collideAndStream() or stream() for the heavy fluid.
    vector<MultiBlockLattice3D<T, DESCRIPTOR>* > blockLattices;
    blockLattices.push_back(&heavyFluid);
    blockLattices.push_back(&lightFluid);
    
    // The argument "constOmegaValues" to the Shan/Chen processor is optional,
    //   and is used for efficiency reasons only. It tells the data processor
    //   that the relaxation times are constant, and that their inverse must be
    //   computed only once.
    std::vector<T> constOmegaValues;
    constOmegaValues.push_back(omega1);
    constOmegaValues.push_back(omega2);
    plint processorLevel = 1;
    integrateProcessingFunctional (
            new ShanChenMultiComponentProcessor3D<T,DESCRIPTOR>(G,constOmegaValues),
            Box3D(0,nx-1,1,ny-2,0,nz-1),
            blockLattices, processorLevel );

    rayleighTaylorSetup(heavyFluid, lightFluid, rho0, rho1, force);
	
    pcout << "Starting simulation" << endl;
    // Main loop over time iterations.
    for (plint iT=0; iT<maxIter; ++iT) {
        if (iT%saveIter==0) {
            writePpms(heavyFluid, lightFluid, iT);
        }

        // Time iteration for the light fluid.
        lightFluid.collideAndStream();
        // Time iteration for the heavy fluid must come after the light fluid,
        //   because the coupling is executed here. You should understand this as follows.
        //   The effect of the coupling is to compute the interaction force between
        //   species, and to precompute density and momentum for each species. This must
        //   be executed *before* collide-and-streaming the fluids, because the collision
        //   step needs to access all these values. In the present case, it is done after
        //   both collide-and-stream step, which means, before the collide-and-stream of
        //   the next iteration (it's the same if you are before or after; the important
        //   point is not to be between the two collide-and-streams of the light and heavy
        //   fluid. As for the initial condition, the coupling is initially performed once
        //   during the function call to heavyFluid.initialize().
        heavyFluid.collideAndStream();

        if (iT%statIter==0) {
            pcout << "Average energy fluid one = "
                  << getStoredAverageEnergy<T>(heavyFluid);
            pcout << ", average energy fluid two = "
                  << getStoredAverageEnergy<T>(lightFluid) << endl;
        }
    }

    return 0;
}

