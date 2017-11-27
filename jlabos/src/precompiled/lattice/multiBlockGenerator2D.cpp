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
 * A 2D multiblock lattice -- template instantiation.
 */

#ifdef COMPILE_2D

#include "multiBlock/multiBlockGenerator2D.h"
#include "multiBlock/multiBlockGenerator2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {

template std::auto_ptr<MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D> > 
 defaultGenerateMultiBlockLattice2D< FLOAT_T,descriptors::DESCRIPTOR_2D > (
        MultiBlockManagement2D const& management, plint unnamedDummyArg );

template
std::auto_ptr<MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D> > clone<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& originalLattice,
        Box2D const& subDomain, bool crop );

template
std::auto_ptr<MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D> > generateMultiBlockLattice<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlock2D const& originalBlock, Box2D const& intersection,
        bool crop );

template
std::auto_ptr<MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D> > generateIntersectMultiBlockLattice<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlock2D const& originalBlock1,
        MultiBlock2D const& originalBlock2, bool crop );

template
std::auto_ptr<MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D> > generateIntersectMultiBlockLattice<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlock2D const& originalBlock1,
        MultiBlock2D const& originalBlock2,
        Box2D const& intersection, bool crop );

template
std::auto_ptr<MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D> > generateJoinMultiBlockLattice<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlock2D const& originalBlock1,
        MultiBlock2D const& originalBlock2 );

template
std::auto_ptr<MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D> > extend<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& originalBlock, Box2D const& addedBlock );

template
std::auto_ptr<MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D> > except<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& originalBlock,
        Box2D const& exceptedBlock );

template
std::auto_ptr<MultiBlockLattice2D<FLOAT_T, descriptors::DESCRIPTOR_2D> > redistribute<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T, descriptors::DESCRIPTOR_2D> const& originalBlock,
        SparseBlockStructure2D const& newBlockStructure );

template
std::auto_ptr<MultiBlockLattice2D<FLOAT_T, descriptors::DESCRIPTOR_2D> > redistribute<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T, descriptors::DESCRIPTOR_2D> const& originalBlock,
        SparseBlockStructure2D const& newBlockStructure,
        Box2D const& intersection, bool crop );

template
std::auto_ptr<MultiBlockLattice2D<FLOAT_T, descriptors::DESCRIPTOR_2D> > align<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T, descriptors::DESCRIPTOR_2D> const& originalBlock,
        MultiBlock2D const& partnerBlock );

template
std::auto_ptr<MultiBlockLattice2D<FLOAT_T, descriptors::DESCRIPTOR_2D> > reparallelize<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T, descriptors::DESCRIPTOR_2D> const& originalBlock );

}

#endif  // COMPILE_2D
