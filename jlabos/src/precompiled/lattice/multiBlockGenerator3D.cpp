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
 * A 3D multiblock lattice -- template instantiation.
 */

#ifdef COMPILE_3D

#include "multiBlock/multiBlockGenerator3D.h"
#include "multiBlock/multiBlockGenerator3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

template std::auto_ptr<MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D> > 
 defaultGenerateMultiBlockLattice3D< FLOAT_T,descriptors::DESCRIPTOR_3D > (
        MultiBlockManagement3D const& management, plint unnamedDummyArg );

template
std::auto_ptr<MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D> > clone<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& originalLattice,
        Box3D const& subDomain, bool crop );

template
std::auto_ptr<MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D> > generateMultiBlockLattice<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlock3D const& originalBlock, Box3D const& intersection,
        bool crop );

template
std::auto_ptr<MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D> > generateIntersectMultiBlockLattice<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlock3D const& originalBlock1,
        MultiBlock3D const& originalBlock2, bool crop );

template
std::auto_ptr<MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D> > generateIntersectMultiBlockLattice<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlock3D const& originalBlock1,
        MultiBlock3D const& originalBlock2,
        Box3D const& intersection, bool crop );

template
std::auto_ptr<MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D> > generateJoinMultiBlockLattice<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlock3D const& originalBlock1,
        MultiBlock3D const& originalBlock2 );

template
std::auto_ptr<MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D> > extend<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& originalBlock, Box3D const& addedBlock );

template
std::auto_ptr<MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D> > except<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& originalBlock,
        Box3D const& exceptedBlock );

template
std::auto_ptr<MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D> > redistribute<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D> const& originalBlock,
        SparseBlockStructure3D const& newBlockStructure );

template
std::auto_ptr<MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D> > redistribute<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D> const& originalBlock,
        SparseBlockStructure3D const& newBlockStructure,
        Box3D const& intersection, bool crop );

template
std::auto_ptr<MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D> > align<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D> const& originalBlock,
        MultiBlock3D const& partnerBlock );

template
std::auto_ptr<MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D> > reparallelize<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D> const& originalBlock );

}

#endif  // COMPILE_3D
