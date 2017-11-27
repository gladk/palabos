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

#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.hh"

namespace plb {

/* *************** Boxed Data Processor functionals ****************** */

template class ReductiveBoxProcessingFunctional2D_N<PRECOMP_T>;
template class MaskedReductiveBoxProcessingFunctional2D_N<PRECOMP_T>;
template class ReductiveBoxProcessingFunctional2D_NN<PRECOMP_T,PRECOMP_T>;
template class MaskedReductiveBoxProcessingFunctional2D_NN<PRECOMP_T,PRECOMP_T>;
template class ReductiveNTensorFieldBoxProcessingFunctional2D<PRECOMP_T>;

template class ReductiveBoxProcessingFunctional2D_S<PRECOMP_T>;
template class ReductiveBoxProcessingFunctional2D_SS<PRECOMP_T,PRECOMP_T>;
template class ReductiveScalarFieldBoxProcessingFunctional2D<PRECOMP_T>;

/* *************** Bounded Boxed Data Processor functionals ****************** */

template class BoundedReductiveBoxProcessingFunctional2D_N<PRECOMP_T>;
template class BoundedMaskedReductiveBoxProcessingFunctional2D_N<PRECOMP_T>;
template class BoundedReductiveBoxProcessingFunctional2D_NN<PRECOMP_T,PRECOMP_T>;
template class BoundedMaskedReductiveBoxProcessingFunctional2D_NN<PRECOMP_T,PRECOMP_T>;
template class BoundedReductiveNTensorFieldBoxProcessingFunctional2D<PRECOMP_T>;

}  // namespace plb
