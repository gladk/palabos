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
 * Helper functions for domain initialization -- header file.
 */
#ifndef SWIG_DATA_INITIALIZER_WRAPPER_3D_H
#define SWIG_DATA_INITIALIZER_WRAPPER_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "core/dynamics.h"
#include "plbWrapper/block/dataInitializerFunctional3D.h"

namespace plb {

/// Initialize scalar-field with the same constant value on each cell.
template<typename T>
void setToConstant(MultiNTensorField3D<T>& field, Box3D domain, T value);

/// Initialize scalar-field with the same constant value on each cell.
template<typename T>
void maskedSetToConstant(MultiNTensorField3D<T>& field, MultiNTensorField3D<int>& mask, Box3D domain, T value);

/// Initialize tensor-field with the same constant tensor/vector on each cell.
template<typename T>
void setToConstant( MultiNTensorField3D<T>& field, Box3D domain, T* in_value, int size );

/// Initialize tensor-field with the same constant tensor/vector on each cell.
template<typename T>
void maskedSetToConstant( MultiNTensorField3D<T>& field, MultiNTensorField3D<int>& mask, Box3D domain, T* in_value, int size );

/// Assign the component "index" of its space coordinate to each cell.
template<typename T>
void setToCoordinate(MultiNTensorField3D<T>& field, Box3D domain, plint index);

/// Assign the component "index" of its space coordinate to each cell.
template<typename T>
void maskedSetToCoordinate(MultiNTensorField3D<T>& field, MultiNTensorField3D<int>& mask, Box3D domain, plint index);

/// Assign its space coordinate to each cell.
template<typename T>
void setToCoordinates(MultiNTensorField3D<T>& field,  Box3D domain);

/// Assign its space coordinate to each cell.
template<typename T>
void maskedSetToCoordinates(MultiNTensorField3D<T>& field, MultiNTensorField3D<int>& mask, Box3D domain);

/// Assign scalar-field to one component of a tensor-field.
template<typename T>
void assignComponent(MultiNTensorField3D<T>& tensorField, int whichComponent,
                     MultiNTensorField3D<T>& scalarField, Box3D domain);

/// Assign scalar-field to one component of a tensor-field.
template<typename T>
void maskedAssignComponent(MultiNTensorField3D<T>& tensorField, int whichComponent,
                           MultiNTensorField3D<T>& scalarField, MultiNTensorField3D<int>& mask, Box3D domain);

}  // namespace plb

#endif  // SWIG_DATA_INITIALIZER_WRAPPER_3D_H
