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
 * Flow around a 3D cylinder inside a channel, with the creation of a von
 * Karman vortex street. This example makes use of bounce-back nodes to
 * describe the shape of the cylinder. The outlet is modeled through a
 * Neumann (zero velocity-gradient) condition.
 */

#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>
#include <map>
#include <memory>
#include <algorithm>

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19Descriptor

void cavitySetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                  IncomprFlowParam<T> const& parameters,
                  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    Box3D topLid = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);
    Box3D everythingButTopLid = Box3D(0, nx-1, 0, ny-2, 0, nz-1);

    // All walls implement a Dirichlet velocity condition.
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    T u = std::sqrt((T)2)/(T)2 * parameters.getLatticeU();
    initializeAtEquilibrium(lattice, everythingButTopLid, (T) 1., Array<T,3>((T)0.,(T)0.,(T)0.) );
    initializeAtEquilibrium(lattice, topLid, (T) 1., Array<T,3>(u,(T)0.,u) );
    setBoundaryVelocity(lattice, topLid, Array<T,3>(u,(T)0.,u) );

    lattice.initialize();
}

void printDynamics(std::vector<Dynamics<T,DESCRIPTOR>*> const& dynamics) {
    for (pluint iDyn=0; iDyn<dynamics.size(); ++iDyn) {
        std::vector<int> chain;
        constructIdChain(*dynamics[iDyn], chain);
        pcout << "Structure is: "
              <<  meta::constructIdNameChain<T,DESCRIPTOR>(chain, " >> ") << std::endl;
        pcout << "Omega=" << dynamics[iDyn]->getOmega() << std::endl;
    }
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    std::vector<char> data;
    std::vector<Dynamics<T,DESCRIPTOR>*> old_dynamics;

    old_dynamics.push_back(new BGKdynamics<T,DESCRIPTOR>(3.14159));
    old_dynamics.push_back(new DensityDirichletBoundaryDynamics<T,DESCRIPTOR,1,-1>(new BGKdynamics<T,DESCRIPTOR>(1.4)));

    pcout << "Old dynamics:" << std::endl;
    printDynamics(old_dynamics);
    pcout << std::endl;

    serialize(old_dynamics, data);
    std::vector<Dynamics<T,DESCRIPTOR>*> new_dynamics;
    generateAndUnserializeDynamics(data, new_dynamics);
    pcout << "New dynamics:" << std::endl;
    printDynamics(new_dynamics);

    BGKdynamics<T,DESCRIPTOR> bgkClone(1.);
    pcout << "Before unserialization, BGK clone has omega " << bgkClone.getOmega() << std::endl;
    unserialize(bgkClone, data);
    pcout << "After unserialization, BGK clone has omega " << bgkClone.getOmega() << std::endl;

    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 400.,  // Re
            100,        // N
            1.,        // lx
            1.,        // ly 
            1.         // lz 
    );
    const T logT     = (T)0.02;
    const T maxT     = (T)20.1;

    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
            nx, ny, nz,
            new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    cavitySetup(lattice, parameters, *boundaryCondition);

    std::map<int,std::string> nameOfDynamics;
    ImageWriter<int>("earth").writeScaledGif(
            "dynamics",
            *extractDynamicsChain(lattice,nameOfDynamics,Box3D(nx/2,nx/2,0,ny-1,0,nz-1)), 600,600 );
    //ImageWriter<int>("air").writeScaledGif("dynamics",
    //                                       *extractTopMostDynamics(lattice), 600,600 );
    for (std::map<int,std::string>::const_iterator it = nameOfDynamics.begin();
         it != nameOfDynamics.end();
         ++it)
    { 
        pcout << it->first << " --> " << it->second << std::endl;
    }
    pcout << std::endl;

    std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > newLattice = copyEntireCells(lattice);
    copyEntireCells(*newLattice, lattice, lattice.getBoundingBox());

    // Main loop over time iterations.
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();

        // At this point, the state of the lattice corresponds to the
        //   discrete time iT+1, and the stored averages are upgraded to time iT.
        if (iT%parameters.nStep(logT)==0) {
            pcout << "; av energy ="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho ="
                  << getStoredAverageDensity<T>(lattice) << endl;
        }
    }
    
    delete boundaryCondition;
}
