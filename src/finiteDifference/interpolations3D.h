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

#ifndef INTERPOLATIONS_3D_H
#define INTERPOLATIONS_3D_H

#include "core/globalDefs.h"
#include "core/geometry3D.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/dataField3D.h"
#include <vector>

namespace plb {

/// Helper function: linear interpolation within one cell.
template<typename T>
void linearInterpolationCoefficients (
        AtomicBlock3D const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights );

template<typename T, plint nDim>
Array<T,nDim> linearInterpolateTensorField (
        TensorField3D<T,nDim>& tensorField, Array<T,3> const& position );

template<typename T>
void linearInterpolateNtensorField (
        NTensorField3D<T>& tensorField, Array<T,3> const& position,
        std::vector<T>& result );

template<typename T, plint nDim>
Array<T,nDim> predictorCorrectorTensorField (
        TensorField3D<T,nDim>& tensorField, Array<T,3> const& position, T scaling );

template<typename T>
Array<T,3> predictorCorrectorNTensorField (
        NTensorField3D<T>& tensorField, Array<T,3> const& position, T scaling );

template<typename T>
void predictorCorrectorRhoBarJ (
        NTensorField3D<T>& rhoBarJ, Array<T,3> const& position,
        bool velIsJ, Array<T,3>& j, T& rhoBar );

template<typename T>
T linearInterpolateScalarField (
        ScalarField3D<T>& scalarField, Array<T,3> const& position );

}  // namespace plb

#endif  // INTERPOLATIONS_3D_H

