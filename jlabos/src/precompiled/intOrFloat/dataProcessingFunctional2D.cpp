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

#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/dataProcessingFunctional2D.hh"

namespace plb {


/* *************** Boxed Data Processor functionals ****************** */

template class BoxProcessingFunctional2D_N<PRECOMP_T>;
template class BoxProcessingFunctional2D_S<PRECOMP_T>;
template class MaskedBoxProcessingFunctional2D_N<PRECOMP_T>;
template class BoxProcessingFunctional2D_NN<PRECOMP_T,int>;
template class BoxProcessingFunctional2D_NN<PRECOMP_T,float>;
template class BoxProcessingFunctional2D_NN<PRECOMP_T,double>;
template class BoxProcessingFunctional2D_SS<PRECOMP_T,int>;
template class BoxProcessingFunctional2D_SS<PRECOMP_T,float>;
template class BoxProcessingFunctional2D_SS<PRECOMP_T,double>;
template class MaskedBoxProcessingFunctional2D_NN<PRECOMP_T,float>;
template class MaskedBoxProcessingFunctional2D_NN<PRECOMP_T,double>;
template class MaskedBoxProcessingFunctional2D_NN<PRECOMP_T,int>;
template class NTensorFieldBoxProcessingFunctional2D<PRECOMP_T>;
template class ScalarFieldBoxProcessingFunctional2D<PRECOMP_T>;
template class MaskedNTensorFieldBoxProcessingFunctional2D<PRECOMP_T>;

/* *************** Bounded Boxed Data Processor functionals ****************** */

template class BoundedBoxProcessingFunctional2D_N<PRECOMP_T>;
template class BoundedMaskedBoxProcessingFunctional2D_N<PRECOMP_T>;
template class BoundedBoxProcessingFunctional2D_NN<PRECOMP_T,PRECOMP_T>;
template class BoundedMaskedBoxProcessingFunctional2D_NN<PRECOMP_T,PRECOMP_T>;
template class BoundedNTensorFieldBoxProcessingFunctional2D<PRECOMP_T>;

}  // namespace plb
