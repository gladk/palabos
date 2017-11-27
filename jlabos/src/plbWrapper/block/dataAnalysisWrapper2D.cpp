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

#include "plbWrapper/block/dataAnalysisWrapper2D.h"
#include "plbWrapper/block/dataAnalysisWrapper2D.hh"
#include "dataProcessors/ntensorAnalysisWrapper2D.h"
#include "dataProcessors/ntensorAnalysisWrapper2D.hh"

namespace plb {

template
void computeAverage (
        MultiNTensorField2D<PRECOMP_T>& vectorField, Box2D domain, PRECOMP_T* result, int size );

template
void maskedComputeAverage (
        MultiNTensorField2D<PRECOMP_T>& vectorField, MultiNTensorField2D<int>& mask, Box2D domain, PRECOMP_T* result, int size );

template
void computeMin (
        MultiNTensorField2D<PRECOMP_T>& vectorField, Box2D domain, PRECOMP_T* result, int size );

template
void maskedComputeMin (
        MultiNTensorField2D<PRECOMP_T>& vectorField, MultiNTensorField2D<int>& mask, Box2D domain, PRECOMP_T* result, int size );

template
void computeMax (
        MultiNTensorField2D<PRECOMP_T>& vectorField, Box2D domain, PRECOMP_T* result, int size );

template
void maskedComputeMax (
        MultiNTensorField2D<PRECOMP_T>& vectorField, MultiNTensorField2D<int>& mask, Box2D domain, PRECOMP_T* result, int size );

template
void computeBoundedAverage (
        MultiNTensorField2D<PRECOMP_T>& vectorField, Box2D domain, PRECOMP_T* result, int size );

template
void maskedComputeBoundedAverage (
        MultiNTensorField2D<PRECOMP_T>& vectorField, MultiNTensorField2D<int>& mask, Box2D domain, PRECOMP_T* result, int size );

template
void copy( MultiNTensorField2D<PRECOMP_T>& field,
           MultiNTensorField2D<int>& convertedField, Box2D domain);

template
void maskedCopy( MultiNTensorField2D<PRECOMP_T>& field,
                 MultiNTensorField2D<int>& convertedField, MultiNTensorField2D<int>& mask,
                 Box2D domain);

template
MultiNTensorField2D<int>* copyConvert (
        MultiNTensorField2D<PRECOMP_T>& field, Box2D domain );

template
MultiNTensorField2D<int>* maskedCopyConvert (
        MultiNTensorField2D<PRECOMP_T>& field, MultiNTensorField2D<int>& mask, Box2D domain );



template
void copy( MultiNTensorField2D<PRECOMP_T>& field,
           MultiNTensorField2D<float>& convertedField, Box2D domain);

template
void maskedCopy( MultiNTensorField2D<PRECOMP_T>& field,
                 MultiNTensorField2D<float>& convertedField, MultiNTensorField2D<int>& mask,
                 Box2D domain);

template
MultiNTensorField2D<float>* copyConvert (
        MultiNTensorField2D<PRECOMP_T>& field, Box2D domain );

template
MultiNTensorField2D<float>* maskedCopyConvert (
        MultiNTensorField2D<PRECOMP_T>& field, MultiNTensorField2D<int>& mask, Box2D domain );



template
void copy( MultiNTensorField2D<PRECOMP_T>& field,
           MultiNTensorField2D<double>& convertedField, Box2D domain);
template
void maskedCopy( MultiNTensorField2D<PRECOMP_T>& field,
                 MultiNTensorField2D<double>& convertedField, MultiNTensorField2D<int>& mask,
                 Box2D domain);

template
MultiNTensorField2D<double>* copyConvert (
        MultiNTensorField2D<PRECOMP_T>& field, Box2D domain );

template
MultiNTensorField2D<double>* maskedCopyConvert (
        MultiNTensorField2D<PRECOMP_T>& field, MultiNTensorField2D<int>& mask, Box2D domain );


template
void extractComponent (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        MultiNTensorField2D<PRECOMP_T>& component,
        Box2D domain, int iComponent );

template
void maskedExtractComponent (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        MultiNTensorField2D<PRECOMP_T>& component, MultiNTensorField2D<int>& mask,
        Box2D domain, int iComponent );

template
MultiNTensorField2D<PRECOMP_T>* extractComponent (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        Box2D domain, int iComponent );

template
MultiNTensorField2D<PRECOMP_T>* maskedExtractComponent (
        MultiNTensorField2D<PRECOMP_T>& tensorField, MultiNTensorField2D<int>& mask,
        Box2D domain, int iComponent );

template
void computeNorm (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        MultiNTensorField2D<PRECOMP_T>& norm, Box2D domain );

template
void maskedComputeNorm (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        MultiNTensorField2D<PRECOMP_T>& norm, MultiNTensorField2D<int>& mask, Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* computeNorm (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedComputeNorm (
        MultiNTensorField2D<PRECOMP_T>& tensorField, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void computeNormSqr (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        MultiNTensorField2D<PRECOMP_T>& normSqr, Box2D domain );

template
void maskedComputeNormSqr (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        MultiNTensorField2D<PRECOMP_T>& normSqr, MultiNTensorField2D<int>& mask, Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* computeNormSqr (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedComputeNormSqr (
        MultiNTensorField2D<PRECOMP_T>& tensorField, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void computeSymmetricTensorNorm (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        MultiNTensorField2D<PRECOMP_T>& norm,
        Box2D domain );

template
void maskedComputeSymmetricTensorNorm (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        MultiNTensorField2D<PRECOMP_T>& norm, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* computeSymmetricTensorNorm (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedComputeSymmetricTensorNorm (
        MultiNTensorField2D<PRECOMP_T>& tensorField, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void computeSymmetricTensorNormSqr (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        MultiNTensorField2D<PRECOMP_T>& normSqr,
        Box2D domain );

template
void maskedComputeSymmetricTensorNormSqr (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        MultiNTensorField2D<PRECOMP_T>& normSqr, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* computeSymmetricTensorNormSqr (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedComputeSymmetricTensorNormSqr (
        MultiNTensorField2D<PRECOMP_T>& tensorField, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void computeSymmetricTensorTrace (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        MultiNTensorField2D<PRECOMP_T>& trace,
        Box2D domain );

template
void maskedComputeSymmetricTensorTrace (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        MultiNTensorField2D<PRECOMP_T>& trace, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* computeSymmetricTensorTrace (
        MultiNTensorField2D<PRECOMP_T>& tensorField,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedComputeSymmetricTensorTrace (
        MultiNTensorField2D<PRECOMP_T>& tensorField, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void computeVorticity (
        MultiNTensorField2D<PRECOMP_T>& velocity,
        MultiNTensorField2D<PRECOMP_T>& vorticity,
        Box2D domain );

template
void maskedComputeVorticity (
        MultiNTensorField2D<PRECOMP_T>& velocity,
        MultiNTensorField2D<PRECOMP_T>& vorticity, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* computeVorticity (
        MultiNTensorField2D<PRECOMP_T>& velocity,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedComputeVorticity (
        MultiNTensorField2D<PRECOMP_T>& velocity, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void computeBulkVorticity (
        MultiNTensorField2D<PRECOMP_T>& velocity,
        MultiNTensorField2D<PRECOMP_T>& vorticity,
        Box2D domain );

template
void maskedComputeBulkVorticity (
        MultiNTensorField2D<PRECOMP_T>& velocity,
        MultiNTensorField2D<PRECOMP_T>& vorticity, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* computeBulkVorticity (
        MultiNTensorField2D<PRECOMP_T>& velocity,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedComputeBulkVorticity (
        MultiNTensorField2D<PRECOMP_T>& velocity, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void computeStrainRate (
        MultiNTensorField2D<PRECOMP_T>& velocity,
        MultiNTensorField2D<PRECOMP_T>& S,
        Box2D domain );

template
void maskedComputeStrainRate (
        MultiNTensorField2D<PRECOMP_T>& velocity,
        MultiNTensorField2D<PRECOMP_T>& S, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* computeStrainRate (
        MultiNTensorField2D<PRECOMP_T>& velocity,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedComputeStrainRate (
        MultiNTensorField2D<PRECOMP_T>& velocity, MultiNTensorField2D<int>& mask,
        Box2D domain );

template void computeBulkStrainRate (
        MultiNTensorField2D<PRECOMP_T>& velocity,
        MultiNTensorField2D<PRECOMP_T>& S,
        Box2D domain );

template
void maskedComputeBulkStrainRate (
        MultiNTensorField2D<PRECOMP_T>& velocity,
        MultiNTensorField2D<PRECOMP_T>& S, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* computeBulkStrainRate (
        MultiNTensorField2D<PRECOMP_T>& velocity,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedComputeBulkStrainRate (
        MultiNTensorField2D<PRECOMP_T>& velocity, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void add (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<PRECOMP_T>& result,
        Box2D domain );

template
void maskedAdd (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<PRECOMP_T>& result, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* add (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedAdd (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void add (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<PRECOMP_T>& result,
        Box2D domain);

template
void maskedAdd (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<PRECOMP_T>& result, MultiNTensorField2D<int>& mask,
        Box2D domain);

template
MultiNTensorField2D<PRECOMP_T>* add (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedAdd (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
void subtract (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<PRECOMP_T>& result,
        Box2D domain );

template
void maskedSubtract (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<PRECOMP_T>& result, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
MultiNTensorField2D<PRECOMP_T>* subtract (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedSubtract (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
void subtract (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<PRECOMP_T>& result,
        Box2D domain );

template
void maskedSubtract (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<PRECOMP_T>& result, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* subtract (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedSubtract (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void subtract (
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& result,
        Box2D domain );

template
void maskedSubtract (
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& result, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* subtract (
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<PRECOMP_T>& A,
        Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedSubtract (
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<PRECOMP_T>& A, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
void multiply (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<PRECOMP_T>& result,
        Box2D domain );

template
void maskedMultiply (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<PRECOMP_T>& result, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
MultiNTensorField2D<PRECOMP_T>* multiply(MultiNTensorField2D<PRECOMP_T>& A,
                                         MultiNTensorField2D<PRECOMP_T>& B, Box2D domain);

template
MultiNTensorField2D<PRECOMP_T>* maskedMultiply(MultiNTensorField2D<PRECOMP_T>& A,
                                               MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask, Box2D domain);


template
void multiply(MultiNTensorField2D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
              MultiNTensorField2D<PRECOMP_T>& result, Box2D domain);

template
void maskedMultiply(MultiNTensorField2D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
              MultiNTensorField2D<PRECOMP_T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template
MultiNTensorField2D<PRECOMP_T>* multiply(MultiNTensorField2D<PRECOMP_T>& A, PRECOMP_T* alpha, int size, Box2D domain);

template
MultiNTensorField2D<PRECOMP_T>* maskedMultiply(MultiNTensorField2D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
                                               MultiNTensorField2D<int>& mask, Box2D domain);

template
void divide (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<PRECOMP_T>& result,
        Box2D domain );

template
void maskedDivide (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<PRECOMP_T>& result, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
MultiNTensorField2D<PRECOMP_T>* divide(MultiNTensorField2D<PRECOMP_T>& A,
                                       MultiNTensorField2D<PRECOMP_T>& B, Box2D domain);

template
MultiNTensorField2D<PRECOMP_T>* maskedDivide(MultiNTensorField2D<PRECOMP_T>& A,
                                             MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask, Box2D domain);


template
void divide(MultiNTensorField2D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
            MultiNTensorField2D<PRECOMP_T>& result, Box2D domain);

template
void maskedDivide(MultiNTensorField2D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
            MultiNTensorField2D<PRECOMP_T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template
MultiNTensorField2D<PRECOMP_T>* divide(MultiNTensorField2D<PRECOMP_T>& A, PRECOMP_T* alpha, int size, Box2D domain);

template
MultiNTensorField2D<PRECOMP_T>* maskedDivide(MultiNTensorField2D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
                                       MultiNTensorField2D<int>& mask, Box2D domain);


template
void divide(PRECOMP_T* alpha, int size, MultiNTensorField2D<PRECOMP_T>& A,
            MultiNTensorField2D<PRECOMP_T>& result, Box2D domain);

template
void maskedDivide(PRECOMP_T* alpha, int size, MultiNTensorField2D<PRECOMP_T>& A,
            MultiNTensorField2D<PRECOMP_T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template
MultiNTensorField2D<PRECOMP_T>* divide (
      PRECOMP_T* alpha, int size, MultiNTensorField2D<PRECOMP_T>& A, Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedDivide (
      PRECOMP_T* alpha, int size, MultiNTensorField2D<PRECOMP_T>& A, MultiNTensorField2D<int>& mask, Box2D domain );



template
void toThePower (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<PRECOMP_T>& result,
        Box2D domain );

template
void maskedToThePower (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<PRECOMP_T>& result, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
MultiNTensorField2D<PRECOMP_T>* toThePower (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedToThePower (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask, Box2D domain );


template
void toThePower(MultiNTensorField2D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
                MultiNTensorField2D<PRECOMP_T>& result, Box2D domain);

template
void maskedToThePower(MultiNTensorField2D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
                MultiNTensorField2D<PRECOMP_T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template
MultiNTensorField2D<PRECOMP_T>* toThePower (
        MultiNTensorField2D<PRECOMP_T>& A, PRECOMP_T* alpha, int size, Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedToThePower (
        MultiNTensorField2D<PRECOMP_T>& A, PRECOMP_T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain );


template
void toThePower(PRECOMP_T* alpha, int size, MultiNTensorField2D<PRECOMP_T>& A,
                MultiNTensorField2D<PRECOMP_T>& result, Box2D domain);

template
void maskedToThePower(PRECOMP_T* alpha, int size, MultiNTensorField2D<PRECOMP_T>& A,
                MultiNTensorField2D<PRECOMP_T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template
MultiNTensorField2D<PRECOMP_T>* toThePower (
      PRECOMP_T* alpha, int size, MultiNTensorField2D<PRECOMP_T>& A, Box2D domain );

template
MultiNTensorField2D<PRECOMP_T>* maskedToThePower (
      PRECOMP_T* alpha, int size, MultiNTensorField2D<PRECOMP_T>& A, MultiNTensorField2D<int>& mask, Box2D domain );


template
void equals (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<int>& result,
        Box2D domain );

template
void maskedEquals (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<int>* equals (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
MultiNTensorField2D<int>* maskedEquals (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void equals (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<int>& result,
        Box2D domain);

template
void maskedEquals (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask,
        Box2D domain);

template
MultiNTensorField2D<int>* equals (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box2D domain );

template
MultiNTensorField2D<int>* maskedEquals (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
void lessThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<int>& result,
        Box2D domain );

template
void maskedLessThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<int>* lessThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
MultiNTensorField2D<int>* maskedLessThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void lessThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<int>& result,
        Box2D domain);

template
void maskedLessThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask,
        Box2D domain);

template
MultiNTensorField2D<int>* lessThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box2D domain );

template
MultiNTensorField2D<int>* maskedLessThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
void lessEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<int>& result,
        Box2D domain );

template
void maskedLessEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<int>* lessEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
MultiNTensorField2D<int>* maskedLessEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void lessEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<int>& result,
        Box2D domain);

template
void maskedLessEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask,
        Box2D domain);

template
MultiNTensorField2D<int>* lessEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box2D domain );

template
MultiNTensorField2D<int>* maskedLessEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
void greaterThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<int>& result,
        Box2D domain );

template
void maskedGreaterThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<int>* greaterThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
MultiNTensorField2D<int>* maskedGreaterThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void greaterThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<int>& result,
        Box2D domain);

template
void maskedGreaterThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask,
        Box2D domain);

template
MultiNTensorField2D<int>* greaterThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box2D domain );

template
MultiNTensorField2D<int>* maskedGreaterThan (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
void greaterEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<int>& result,
        Box2D domain );

template
void maskedGreaterEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<int>* greaterEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
MultiNTensorField2D<int>* maskedGreaterEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void greaterEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<int>& result,
        Box2D domain);

template
void maskedGreaterEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask,
        Box2D domain);

template
MultiNTensorField2D<int>* greaterEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box2D domain );

template
MultiNTensorField2D<int>* maskedGreaterEqual (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
void logicalAnd (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<int>& result,
        Box2D domain );

template
void maskedLogicalAnd (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<int>* logicalAnd (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
MultiNTensorField2D<int>* maskedLogicalAnd (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void logicalOr (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<int>& result,
        Box2D domain );

template
void maskedLogicalOr (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
MultiNTensorField2D<int>* logicalOr (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
MultiNTensorField2D<int>* maskedLogicalOr (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void negate (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<int>& result,
        Box2D domain);

template
void maskedNegate (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask,
        Box2D domain);

template
MultiNTensorField2D<int>* negate (
        MultiNTensorField2D<PRECOMP_T>& A,
        Box2D domain );

template
MultiNTensorField2D<int>* maskedNegate (
        MultiNTensorField2D<PRECOMP_T>& A, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
void assignInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
void maskedAssignInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void assignInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box2D domain );

template
void maskedAssignInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void addInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
void maskedAddInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void addInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box2D domain );

template
void maskedAddInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void subtractInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
void maskedSubtractInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void subtractInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box2D domain );

template
void maskedSubtractInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void multiplyInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
void maskedMultiplyInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void multiplyInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box2D domain );

template
void maskedMultiplyInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void divideInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
void maskedDivideInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void divideInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box2D domain );

template
void maskedDivideInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
void toThePowerInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B,
        Box2D domain );

template
void maskedToThePowerInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        MultiNTensorField2D<PRECOMP_T>& B, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void toThePowerInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box2D domain );

template
void maskedToThePowerInPlace (
        MultiNTensorField2D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField2D<int>& mask,
        Box2D domain );

}  // namespace plb
