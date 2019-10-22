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
 * Helper functions for domain initialization -- generic implementation.
 */
#ifndef SWIG_DATA_INITIALIZER_WRAPPER_3D_HH
#define SWIG_DATA_INITIALIZER_WRAPPER_3D_HH

#include "plbWrapper/block/dataInitializerWrapper3D.h"
#include "core/cell.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T>
void setToConstant(MultiNTensorField3D<T>& field, Box3D domain, T value) {
    std::vector<T> wrapValue(1);
    wrapValue[0] = value;
    applyProcessingFunctional(new IniConstNTensorFunctional3D<T>(wrapValue), domain, field);
}

template<typename T>
void maskedSetToConstant(MultiNTensorField3D<T>& field, MultiNTensorField3D<int>& mask, Box3D domain, T value) {
    std::vector<T> wrapValue(1);
    wrapValue[0] = value;
    applyProcessingFunctional(new MaskedIniConstNTensorFunctional3D<T>(wrapValue), domain, field, mask);
}

template<typename T>
void setToConstant( MultiNTensorField3D<T>& field, Box3D domain,
                    T* in_value, int size )
{
    std::vector<T> wrapValue(in_value, in_value+size);
    applyProcessingFunctional (
            new IniConstNTensorFunctional3D<T>(wrapValue), domain, field );
}

template<typename T>
void maskedSetToConstant( MultiNTensorField3D<T>& field, MultiNTensorField3D<int>& mask, Box3D domain,
                          T* in_value, int size )
{
    std::vector<T> wrapValue(in_value, in_value+size);
    applyProcessingFunctional (
            new MaskedIniConstNTensorFunctional3D<T>(wrapValue), domain, field, mask );
}

template<typename T>
void setToCoordinate(MultiNTensorField3D<T>& field, Box3D domain, plint index) {
    applyProcessingFunctional(new SetToNCoordinateFunctional3D<T>(index), domain, field);
}

template<typename T>
void maskedSetToCoordinate(MultiNTensorField3D<T>& field, MultiNTensorField3D<int>& mask, Box3D domain, plint index) {
    applyProcessingFunctional(new MaskedSetToNCoordinateFunctional3D<T>(index), domain, field, mask);
}

template<typename T>
void setToCoordinates(MultiNTensorField3D<T>& field, Box3D domain) {
    applyProcessingFunctional(new SetToNCoordinatesFunctional3D<T>, domain, field);
}

template<typename T>
void maskedSetToCoordinates(MultiNTensorField3D<T>& field, MultiNTensorField3D<int>& mask, Box3D domain) {
    applyProcessingFunctional(new MaskedSetToNCoordinatesFunctional3D<T>, domain, field, mask);
}

template<typename T>
void assignComponent(MultiNTensorField3D<T>& tensorField, int whichComponent,
                     MultiNTensorField3D<T>& scalarField, Box3D domain)
{
    applyProcessingFunctional(new SetNTensorComponentFunctional3D<T>(whichComponent),
                              domain, scalarField, tensorField);
}

template<typename T>
void maskedAssignComponent(MultiNTensorField3D<T>& tensorField, int whichComponent,
                           MultiNTensorField3D<T>& scalarField, MultiNTensorField3D<int>& mask, Box3D domain)
{
    applyProcessingFunctional(new MaskedSetNTensorComponentFunctional3D<T>(whichComponent),
                              domain, scalarField, tensorField, mask);
}

}  // namespace plb

#endif  // SWIG_DATA_INITIALIZER_WRAPPER_3D_HH
