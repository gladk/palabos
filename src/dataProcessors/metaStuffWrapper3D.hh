/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2013 FlowKit Sarl
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
 * Helper functions for domain initialization -- header file.
 */
#ifndef META_STUFF_WRAPPER_3D_HH
#define META_STUFF_WRAPPER_3D_HH

#include "dataProcessors/metaStuffWrapper3D.h"
#include "dataProcessors/metaStuffFunctional3D.h"
#include "dataProcessors/dataInitializerWrapper3D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "multiBlock/multiContainerBlock3D.h"
#include "core/dynamicsIdentifiers.h"
#include <set>


namespace plb {

template<typename T, template<typename U> class Descriptor>
void extractTopMostDynamics(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<int>& dynamicsId,
                            Box3D domain )
{
    applyProcessingFunctional (
            new ExtractTopMostDynamicsFunctional3D<T,Descriptor>(),
            domain, lattice, dynamicsId );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr< MultiScalarField3D<int> > extractTopMostDynamics (
                                             MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain )
{
    MultiScalarField3D<int>* dynamicsId = new MultiScalarField3D<int>(lattice, domain);
    extractTopMostDynamics(lattice, *dynamicsId, domain);
    return std::auto_ptr<MultiScalarField3D<int> >(dynamicsId);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr< MultiScalarField3D<int> > extractTopMostDynamics (
                                             MultiBlockLattice3D<T,Descriptor>& lattice )
{
    return extractTopMostDynamics(lattice, lattice.getBoundingBox());
}


template<typename T, template<typename U> class Descriptor>
void extractBottomMostDynamics(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<int>& dynamicsId,
                            Box3D domain )
{
    applyProcessingFunctional (
            new ExtractBottomMostDynamicsFunctional3D<T,Descriptor>(),
            domain, lattice, dynamicsId );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr< MultiScalarField3D<int> > extractBottomMostDynamics (
                                             MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain )
{
    MultiScalarField3D<int>* dynamicsId = new MultiScalarField3D<int>(lattice, domain);
    extractBottomMostDynamics(lattice, *dynamicsId, domain);
    return std::auto_ptr<MultiScalarField3D<int> >(dynamicsId);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr< MultiScalarField3D<int> > extractBottomMostDynamics (
                                             MultiBlockLattice3D<T,Descriptor>& lattice )
{
    return extractBottomMostDynamics(lattice, lattice.getBoundingBox());
}


template<typename T, template<typename U> class Descriptor>
void uniqueDynamicsChains (
        MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
        std::vector<std::vector<int> >& chains, pluint& maxChainLength )
{
    MultiContainerBlock3D container(lattice);
    std::vector<MultiBlock3D*> latticeAndContainer, containerVector;
    latticeAndContainer.push_back(&lattice);
    latticeAndContainer.push_back(&container);
    containerVector.push_back(&container);

    StoreDynamicsFunctional3D<T,Descriptor> SDFunctional;
    applyProcessingFunctional(SDFunctional, domain, latticeAndContainer);
    maxChainLength = SDFunctional.getMaxChainLength();

    std::vector<int> nextMaximum(maxChainLength);
    for (pluint i=0; i<nextMaximum.size(); ++i) {
        nextMaximum[i] = -1;
    }
    std::vector<int> emptyVector(nextMaximum);
    std::vector<std::vector<int> > maxima;
    do {
        IterateDynamicsFunctional3D functional(nextMaximum);
        applyProcessingFunctional(functional, domain, containerVector);
        nextMaximum = functional.getNextMaximum();
        if (!vectorEquals(nextMaximum, emptyVector)) {
            maxima.push_back(nextMaximum);
        }
    }
    while (!vectorEquals(nextMaximum, emptyVector));
    chains.swap(maxima);
}

template<typename T, template<typename U> class Descriptor>
void uniqueDynamicsIds (
        MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain, std::vector<int>& ids )
{
    std::vector<std::vector<int> > chains;
    pluint maxChainLength;
    uniqueDynamicsChains(lattice, domain, chains, maxChainLength);
    std::set<int> idSet;
    for (pluint iChain=0; iChain<chains.size(); ++iChain) {
        idSet.insert(chains[iChain].begin(), chains[iChain].end());
    }
    idSet.erase(-1);
    std::vector<int>(idSet.begin(), idSet.end()).swap(ids);
}

template<typename T, template<typename U> class Descriptor>
void extractDynamicsChain(MultiBlockLattice3D<T,Descriptor>& lattice,
                          MultiScalarField3D<int>& dynamicsId,
                          std::map<int,std::string>& nameOfDynamics, Box3D domain )
{
    std::vector<std::vector<int> > chains;
    pluint maxChainLength;
    uniqueDynamicsChains(lattice, domain, chains, maxChainLength);
    nameOfDynamics.clear();
    typename ExtractDynamicsChainFunctional3D<T,Descriptor>::DMap dynamicsMap;
    for (pluint iChain=0; iChain<chains.size(); ++iChain) {
        dynamicsMap[chains[iChain]] = iChain;
        nameOfDynamics[iChain] = meta::constructIdNameChain<T,Descriptor>(chains[iChain], " >> ");
    }
    setToConstant(dynamicsId, dynamicsId.getBoundingBox(), -1);
    applyProcessingFunctional (
            new ExtractDynamicsChainFunctional3D<T,Descriptor>(dynamicsMap, maxChainLength),
            domain, lattice, dynamicsId );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr< MultiScalarField3D<int> > extractDynamicsChain (
            MultiBlockLattice3D<T,Descriptor>& lattice,
            std::map<int,std::string>& nameOfDynamics, Box3D domain )
{
    MultiScalarField3D<int>* dynamicsId = new MultiScalarField3D<int>(lattice, domain);
    extractDynamicsChain(lattice, *dynamicsId, nameOfDynamics, domain);
    return std::auto_ptr<MultiScalarField3D<int> >(dynamicsId);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr< MultiScalarField3D<int> > extractDynamicsChain (
            MultiBlockLattice3D<T,Descriptor>& lattice,
            std::map<int,std::string>& nameOfDynamics )
{
    return extractDynamicsChain(lattice, nameOfDynamics, lattice.getBoundingBox());
}


template<typename T, template<typename U> class Descriptor>
void copyEntireCells( MultiBlockLattice3D<T,Descriptor>& sourceLattice,
                      MultiBlockLattice3D<T,Descriptor>& destinationLattice,
                      Box3D domain )
{
    applyProcessingFunctional(new AssignEntireCellFunctional3D<T,Descriptor>, domain,
                              sourceLattice, destinationLattice);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr< MultiBlockLattice3D<T,Descriptor> > copyEntireCells (
            MultiBlockLattice3D<T,Descriptor>& lattice )
{
    return copyEntireCells(lattice, lattice.getBoundingBox());
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr< MultiBlockLattice3D<T,Descriptor> > copyEntireCells (
            MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain )
{
    MultiBlockLattice3D<T,Descriptor>* newLattice = new MultiBlockLattice3D<T,Descriptor>(lattice, domain);
    copyEntireCells(lattice, *newLattice, domain);
    return std::auto_ptr<MultiBlockLattice3D<T,Descriptor> >(newLattice);
}

}  // namespace plb

#endif  // META_STUFF_WRAPPER_3D_HH
