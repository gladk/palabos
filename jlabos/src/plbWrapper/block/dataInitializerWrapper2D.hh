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
#ifndef SWIG_DATA_INITIALIZER_WRAPPER_2D_HH
#define SWIG_DATA_INITIALIZER_WRAPPER_2D_HH

#include "plbWrapper/block/dataInitializerWrapper2D.h"
#include "core/cell.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T>
void setToConstant(MultiNTensorField2D<T>& field, Box2D domain, T value) {
    std::vector<T> wrapValue(1);
    wrapValue[0] = value;
    applyProcessingFunctional(new IniConstNTensorFunctional2D<T>(wrapValue), domain, field);
}

template<typename T>
void maskedSetToConstant(MultiNTensorField2D<T>& field, MultiNTensorField2D<int>& mask, Box2D domain, T value) {
    std::vector<T> wrapValue(1);
    wrapValue[0] = value;
    applyProcessingFunctional(new MaskedIniConstNTensorFunctional2D<T>(wrapValue), domain, field, mask);
}

template<typename T>
void setToConstant( MultiNTensorField2D<T>& field, Box2D domain,
                    T* in_value, int size )
{
    std::vector<T> wrapValue(in_value, in_value+size);
    applyProcessingFunctional (
            new IniConstNTensorFunctional2D<T>(wrapValue), domain, field );
}

template<typename T>
void maskedSetToConstant( MultiNTensorField2D<T>& field, MultiNTensorField2D<int>& mask, Box2D domain,
                          T* in_value, int size )
{
    std::vector<T> wrapValue(in_value, in_value+size);
    applyProcessingFunctional (
            new MaskedIniConstNTensorFunctional2D<T>(wrapValue), domain, field, mask );
}

template<typename T>
void setToCoordinate(MultiNTensorField2D<T>& field, Box2D domain, plint index) {
    applyProcessingFunctional(new SetToNCoordinateFunctional2D<T>(index), domain, field);
}

template<typename T>
void maskedSetToCoordinate(MultiNTensorField2D<T>& field, MultiNTensorField2D<int>& mask, Box2D domain, plint index) {
    applyProcessingFunctional(new MaskedSetToNCoordinateFunctional2D<T>(index), domain, field, mask);
}

template<typename T>
void setToCoordinates(MultiNTensorField2D<T>& field, Box2D domain) {
    applyProcessingFunctional(new SetToNCoordinatesFunctional2D<T>, domain, field);
}

template<typename T>
void maskedSetToCoordinates(MultiNTensorField2D<T>& field, MultiNTensorField2D<int>& mask, Box2D domain) {
    applyProcessingFunctional(new MaskedSetToNCoordinatesFunctional2D<T>, domain, field, mask);
}

template<typename T>
void assignComponent(MultiNTensorField2D<T>& tensorField, int whichComponent,
                     MultiNTensorField2D<T>& scalarField, Box2D domain)
{
    applyProcessingFunctional(new SetNTensorComponentFunctional2D<T>(whichComponent),
                              domain, scalarField, tensorField);
}

template<typename T>
void maskedAssignComponent(MultiNTensorField2D<T>& tensorField, int whichComponent,
                           MultiNTensorField2D<T>& scalarField, MultiNTensorField2D<int>& mask, Box2D domain)
{
    applyProcessingFunctional(new MaskedSetNTensorComponentFunctional2D<T>(whichComponent),
                              domain, scalarField, tensorField, mask);
}

}  // namespace plb

#endif  // SWIG_DATA_INITIALIZER_WRAPPER_2D_HH
