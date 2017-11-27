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

#ifdef COMPILE_2D

#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {

/* *************** Boxed Data Processor functionals ****************** */

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        ReductiveBoxProcessingFunctional2D_L<FLOAT_T,descriptors::DESCRIPTOR_2D>& functional,
        Box2D domain, MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_2D,
                               FLOAT_T, descriptors::DESCRIPTOR_2D> (
        ReductiveBoxProcessingFunctional2D_LL<FLOAT_T,descriptors::DESCRIPTOR_2D,
                                              FLOAT_T,descriptors::DESCRIPTOR_2D>& functional,
        Box2D domain,
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice1,
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice2 );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_2D, FLOAT_T> (
        ReductiveBoxProcessingFunctional2D_LN< FLOAT_T,descriptors::DESCRIPTOR_2D,
                                               FLOAT_T >& functional,
        Box2D domain,
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& field );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_2D, FLOAT_T> (
        ReductiveBoxProcessingFunctional2D_LS< FLOAT_T,descriptors::DESCRIPTOR_2D,
                                               FLOAT_T >& functional,
        Box2D domain,
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiScalarField2D<FLOAT_T>& field );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_2D, FLOAT_T> (
        MaskedReductiveBoxProcessingFunctional2D_LN< FLOAT_T,descriptors::DESCRIPTOR_2D,
                                                     FLOAT_T >& functional,
        Box2D domain,
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& field,
        MultiNTensorField2D<int>& mask );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        ReductiveLatticeBoxProcessingFunctional2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& functional,
        Box2D domain, std::vector<MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>*> lattices );


/* *************** Bounded Boxed Data Processor functionals ****************** */

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        BoundedReductiveBoxProcessingFunctional2D_L<FLOAT_T,descriptors::DESCRIPTOR_2D>& functional,
        Box2D domain, MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice, plint boundaryWidth );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_2D,
                               FLOAT_T, descriptors::DESCRIPTOR_2D> (
        BoundedReductiveBoxProcessingFunctional2D_LL <
            FLOAT_T,descriptors::DESCRIPTOR_2D,
            FLOAT_T,descriptors::DESCRIPTOR_2D >& functional,
        Box2D domain,
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice1,
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice2,
        plint boundaryWidth );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_2D, FLOAT_T> (
        BoundedReductiveBoxProcessingFunctional2D_LN <
            FLOAT_T, descriptors::DESCRIPTOR_2D, FLOAT_T >& functional,
        Box2D domain,
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& field,
        plint boundaryWidth );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_2D, FLOAT_T> (
        BoundedMaskedReductiveBoxProcessingFunctional2D_LN <
            FLOAT_T, descriptors::DESCRIPTOR_2D, FLOAT_T >& functional,
        Box2D domain,
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& field,
        MultiNTensorField2D<int>& mask,
        plint boundaryWidth );

template
void applyProcessingFunctional<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        BoundedReductiveLatticeBoxProcessingFunctional2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& functional,
        Box2D domain, std::vector<MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>*> lattices,
        plint boundaryWidth);

}  // namespace plb

#endif  // COMPILE_2D
