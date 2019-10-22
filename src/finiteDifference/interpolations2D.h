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

#ifndef INTERPOLATIONS_2D_H
#define INTERPOLATIONS_2D_H

#include "core/globalDefs.h"
#include "core/geometry2D.h"
#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/dataField2D.h"
#include <vector>

namespace plb {

/// Helper function: linear interpolation within one cell.
template<typename T>
void linearInterpolationCoefficients (
        AtomicBlock2D const& block, Array<T,2> const& position );

template<typename T, plint nDim>
Array<T,nDim> linearInterpolateTensorField (
        TensorField2D<T,nDim>& tensorField, Array<T,2> const& position );

template<typename T, plint nDim>
Array<T,nDim> predictorCorrectorTensorField (
        TensorField2D<T,nDim>& tensorField, Array<T,2> const& position, T scaling );

template<typename T>
Array<T,2> predictorCorrectorNTensorField (
        NTensorField2D<T>& tensorField, Array<T,2> const& position, T scaling );

        
template<typename T>
void predictorCorrectorRhoBarJ (
        NTensorField2D<T>& rhoBarJ, Array<T,2> const& position,
        bool velIsJ, Array<T,2>& j, T& rhoBar );
        
template<typename T>
T linearInterpolateScalarField (
        ScalarField2D<T>& scalarField, Array<T,2> const& position );

}  // namespace plb

#endif  // INTERPOLATIONS_2D_H
