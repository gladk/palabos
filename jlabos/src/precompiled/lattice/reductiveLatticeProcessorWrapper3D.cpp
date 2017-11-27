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

#ifdef COMPILE_3D

#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

/* *************** Boxed Data Processor functionals ****************** */

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        ReductiveBoxProcessingFunctional3D_L<FLOAT_T,descriptors::DESCRIPTOR_3D>& functional,
        Box3D domain, BlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_3D,
                               FLOAT_T, descriptors::DESCRIPTOR_3D> (
        ReductiveBoxProcessingFunctional3D_LL<FLOAT_T,descriptors::DESCRIPTOR_3D,
                                     FLOAT_T,descriptors::DESCRIPTOR_3D>& functional,
        Box3D domain,
        BlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice1,
        BlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice2 );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        ReductiveLatticeBoxProcessingFunctional3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& functional,
        Box3D domain, std::vector<BlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>*> lattices );


/* *************** Dotted Data Processor functionals ***************** */

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        ReductiveDotProcessingFunctional3D_L<FLOAT_T,descriptors::DESCRIPTOR_3D>& functional,
        DotList3D const& dotList, BlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_3D,
                               FLOAT_T, descriptors::DESCRIPTOR_3D> (
        ReductiveDotProcessingFunctional3D_LL<FLOAT_T,descriptors::DESCRIPTOR_3D,
                                     FLOAT_T,descriptors::DESCRIPTOR_3D>& functional,
        DotList3D const& dotList,
        BlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice1,
        BlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice2 );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        ReductiveLatticeDotProcessingFunctional3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& functional,
        DotList3D const& dotList, std::vector<BlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>*> lattices );


/* *************** Bounded Boxed Data Processor functionals ****************** */

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        BoundedReductiveBoxProcessingFunctional3D_L<FLOAT_T,descriptors::DESCRIPTOR_3D>& functional,
        Box3D domain, BlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice, plint boundaryWidth );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_3D,
                               FLOAT_T, descriptors::DESCRIPTOR_3D> (
        BoundedReductiveBoxProcessingFunctional3D_LL <
            FLOAT_T,descriptors::DESCRIPTOR_3D,
            FLOAT_T,descriptors::DESCRIPTOR_3D >& functional,
        Box3D domain,
        BlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice1,
        BlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice2,
        plint boundaryWidth );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        BoundedReductiveLatticeBoxProcessingFunctional3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& functional,
        Box3D domain, std::vector<BlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>*> lattices,
        plint boundaryWidth);

}  // namespace plb

#endif  // COMPILE_3D
