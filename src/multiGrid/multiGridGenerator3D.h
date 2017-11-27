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

/* Main author: Daniel Lagrava
 **/

/** \file
 * Various factories that use a multiGridManagement3D -- Header file
 */

#ifndef MULTI_GRID_GENERATOR_3D_H
#define MULTI_GRID_GENERATOR_3D_H

#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiGrid/multiGridManagement3D.h"

namespace plb {

/// Compute the reduced bulk of the domain for fine to coarse copy
inline Box3D computeCopyReducedBulk(Box3D domain);
inline void computeCopyEdges(Box3D domain, std::vector<Box3D>& edges);
inline void computeCopyCorners(Box3D domain, std::vector<Box3D>& corners);

template<typename T>
void computeFilteringIndicesEdges(Box3D domain, std::vector<std::vector<plint> >& indices);


template<typename T, template<typename U> class Descriptor>
std::vector<MultiBlockLattice3D<T,Descriptor>*> generateLattices(
                MultiGridManagement3D management,
                std::vector<Dynamics<T,Descriptor>*> backgroundDynamics,
                std::vector<BlockCommunicator3D*> communicators,
                std::vector<CombinedStatistics*> combinedStatistics );

template<typename T, template<typename U> class Descriptor>
std::vector<MultiBlockLattice3D<T,Descriptor>*> generateLattices(
                MultiGridManagement3D management,
                std::vector<Dynamics<T,Descriptor>*> backgroundDynamics);


template<typename T, template<typename U> class Descriptor>
void createInterfaces(std::vector<MultiBlockLattice3D<T,Descriptor>*>& multiBlocks,
                      MultiGridManagement3D management);


template<typename T, template<typename U> class Descriptor>
void createCoarseGridInterface (
        plint coarseLevel, Box3D coarseGridInterface, std::vector<MultiBlockLattice3D<T,Descriptor>*>& multiBlocks );


template<typename T, template<typename U> class Descriptor>
void createFineGridInterface (
        plint coarseLevel, Box3D fineGridInterface, 
        std::vector<MultiBlockLattice3D<T,Descriptor>*>& multiBlocks );
    



template<typename T>
std::vector<MultiScalarField3D<T>*> generateScalarFields(
                        MultiGridManagement3D const& management,
                        std::vector<BlockCommunicator3D*> communicators,
                        std::vector<CombinedStatistics*> combinedStatistics );


template<typename T, int nDim>
std::vector<MultiTensorField3D<T,nDim> *> generateTensorFields (
        MultiGridManagement3D const& management,
        std::vector<BlockCommunicator3D*> communicators,
        std::vector<CombinedStatistics*> combinedStatistics );


} // namespace plb

#endif  // MULTI_GRID_GENERATOR_3D_H

