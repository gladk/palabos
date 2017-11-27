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

#include "atomicBlock/reductiveDataProcessorWrapper2D.h"
#include "atomicBlock/reductiveDataProcessorWrapper2D.hh"

namespace plb {

/* *************** Boxed Data Processor functionals ****************** */

template
void applyProcessingFunctional<PRECOMP_T> (
        ReductiveBoxProcessingFunctional2D_N<PRECOMP_T>& functional,
        Box2D domain, NTensorField2D<PRECOMP_T>& field );

template
void applyProcessingFunctional<PRECOMP_T,PRECOMP_T> (
        ReductiveBoxProcessingFunctional2D_NN<PRECOMP_T,PRECOMP_T>& functional,
        Box2D domain,
        NTensorField2D<PRECOMP_T>& field1,
        NTensorField2D<PRECOMP_T>& field2 );

template
void applyProcessingFunctional<PRECOMP_T> (
        ReductiveNTensorFieldBoxProcessingFunctional2D<PRECOMP_T>& functional,
        Box2D domain, std::vector<NTensorField2D<PRECOMP_T>*> fields );

template
void applyProcessingFunctional<PRECOMP_T> (
        ReductiveBoxProcessingFunctional2D_S<PRECOMP_T>& functional,
        Box2D domain, ScalarField2D<PRECOMP_T>& field );

template
void applyProcessingFunctional<PRECOMP_T,PRECOMP_T> (
        ReductiveBoxProcessingFunctional2D_SS<PRECOMP_T,PRECOMP_T>& functional,
        Box2D domain,
        ScalarField2D<PRECOMP_T>& field1,
        ScalarField2D<PRECOMP_T>& field2 );

template
void applyProcessingFunctional<PRECOMP_T> (
        ReductiveScalarFieldBoxProcessingFunctional2D<PRECOMP_T>& functional,
        Box2D domain, std::vector<ScalarField2D<PRECOMP_T>*> fields );

/* *************** Bounded Boxed Data Processor functionals ****************** */

template
void applyProcessingFunctional<PRECOMP_T> (
        BoundedReductiveBoxProcessingFunctional2D_N<PRECOMP_T>& functional,
        Box2D domain, NTensorField2D<PRECOMP_T>& field, plint boundaryWidth );

template
void applyProcessingFunctional<PRECOMP_T,PRECOMP_T> (
        BoundedReductiveBoxProcessingFunctional2D_NN<PRECOMP_T,PRECOMP_T>& functional,
        Box2D domain,
        NTensorField2D<PRECOMP_T>& field1,
        NTensorField2D<PRECOMP_T>& field2,
        plint boundaryWidth );

}  // namespace plb
