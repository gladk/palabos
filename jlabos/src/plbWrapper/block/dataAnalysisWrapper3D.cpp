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

#include "plbWrapper/block/dataAnalysisWrapper3D.h"
#include "plbWrapper/block/dataAnalysisWrapper3D.hh"
#include "dataProcessors/ntensorAnalysisWrapper3D.h"
#include "dataProcessors/ntensorAnalysisWrapper3D.hh"

namespace plb {

template
void computeAverage (
        MultiNTensorField3D<PRECOMP_T>& vectorField, Box3D domain, PRECOMP_T* result, int size );

template
void maskedComputeAverage (
        MultiNTensorField3D<PRECOMP_T>& vectorField, MultiNTensorField3D<int>& mask,
        Box3D domain, PRECOMP_T* result, int size );

template
void computeMin (
        MultiNTensorField3D<PRECOMP_T>& vectorField, Box3D domain, PRECOMP_T* result, int size );

template
void maskedComputeMin (
        MultiNTensorField3D<PRECOMP_T>& vectorField, MultiNTensorField3D<int>& mask,
        Box3D domain, PRECOMP_T* result, int size );

template
void computeMax (
        MultiNTensorField3D<PRECOMP_T>& vectorField, Box3D domain, PRECOMP_T* result, int size );

template
void maskedComputeMax (
        MultiNTensorField3D<PRECOMP_T>& vectorField, MultiNTensorField3D<int>& mask, Box3D domain, PRECOMP_T* result, int size );

template
void computeBoundedAverage (
        MultiNTensorField3D<PRECOMP_T>& vectorField, Box3D domain, PRECOMP_T* result, int size );

template
void maskedComputeBoundedAverage (
        MultiNTensorField3D<PRECOMP_T>& vectorField, MultiNTensorField3D<int>& mask, Box3D domain, PRECOMP_T* result, int size );

template
void copy( MultiNTensorField3D<PRECOMP_T>& field,
           MultiNTensorField3D<int>& convertedField,
           Box3D domain);

template
void maskedCopy( MultiNTensorField3D<PRECOMP_T>& field,
                 MultiNTensorField3D<int>& convertedField, MultiNTensorField3D<int>& mask,
                 Box3D domain);

template
MultiNTensorField3D<int>* copyConvert (
        MultiNTensorField3D<PRECOMP_T>& field, Box3D domain );

template
MultiNTensorField3D<int>* maskedCopyConvert (
        MultiNTensorField3D<PRECOMP_T>& field, MultiNTensorField3D<int>& mask, Box3D domain );



template
void copy( MultiNTensorField3D<PRECOMP_T>& field,
           MultiNTensorField3D<float>& convertedField, Box3D domain);

template
void maskedCopy( MultiNTensorField3D<PRECOMP_T>& field,
                 MultiNTensorField3D<float>& convertedField, MultiNTensorField3D<int>& mask,
                 Box3D domain);

template
MultiNTensorField3D<float>* copyConvert (
        MultiNTensorField3D<PRECOMP_T>& field, Box3D domain );

template
MultiNTensorField3D<float>* maskedCopyConvert (
        MultiNTensorField3D<PRECOMP_T>& field, MultiNTensorField3D<int>& mask, Box3D domain );



template
void copy( MultiNTensorField3D<PRECOMP_T>& field,
           MultiNTensorField3D<double>& convertedField, Box3D domain);
template
void maskedCopy( MultiNTensorField3D<PRECOMP_T>& field,
                 MultiNTensorField3D<double>& convertedField, MultiNTensorField3D<int>& mask,
                 Box3D domain);

template
MultiNTensorField3D<double>* copyConvert (
        MultiNTensorField3D<PRECOMP_T>& field, Box3D domain );

template
MultiNTensorField3D<double>* maskedCopyConvert (
        MultiNTensorField3D<PRECOMP_T>& field, MultiNTensorField3D<int>& mask, Box3D domain );


template
void extractComponent (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        MultiNTensorField3D<PRECOMP_T>& component,
        Box3D domain, int iComponent );

template
void maskedExtractComponent (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        MultiNTensorField3D<PRECOMP_T>& component, MultiNTensorField3D<int>& mask,
        Box3D domain, int iComponent );

template
MultiNTensorField3D<PRECOMP_T>* extractComponent (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        Box3D domain, int iComponent );

template
MultiNTensorField3D<PRECOMP_T>* maskedExtractComponent (
        MultiNTensorField3D<PRECOMP_T>& tensorField, MultiNTensorField3D<int>& mask,
        Box3D domain, int iComponent );

template
void computeNorm (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        MultiNTensorField3D<PRECOMP_T>& norm, Box3D domain );

template
void maskedComputeNorm (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        MultiNTensorField3D<PRECOMP_T>& norm, MultiNTensorField3D<int>& mask, Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* computeNorm (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedComputeNorm (
        MultiNTensorField3D<PRECOMP_T>& tensorField, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void computeNormSqr (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        MultiNTensorField3D<PRECOMP_T>& normSqr, Box3D domain );

template
void maskedComputeNormSqr (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        MultiNTensorField3D<PRECOMP_T>& normSqr, MultiNTensorField3D<int>& mask, Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* computeNormSqr (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedComputeNormSqr (
        MultiNTensorField3D<PRECOMP_T>& tensorField, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void computeSymmetricTensorNorm (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        MultiNTensorField3D<PRECOMP_T>& norm,
        Box3D domain );

template
void maskedComputeSymmetricTensorNorm (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        MultiNTensorField3D<PRECOMP_T>& norm, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* computeSymmetricTensorNorm (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedComputeSymmetricTensorNorm (
        MultiNTensorField3D<PRECOMP_T>& tensorField, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void computeSymmetricTensorNormSqr (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        MultiNTensorField3D<PRECOMP_T>& normSqr,
        Box3D domain );

template
void maskedComputeSymmetricTensorNormSqr (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        MultiNTensorField3D<PRECOMP_T>& normSqr, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* computeSymmetricTensorNormSqr (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedComputeSymmetricTensorNormSqr (
        MultiNTensorField3D<PRECOMP_T>& tensorField, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void computeSymmetricTensorTrace (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        MultiNTensorField3D<PRECOMP_T>& trace,
        Box3D domain );

template
void maskedComputeSymmetricTensorTrace (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        MultiNTensorField3D<PRECOMP_T>& trace, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* computeSymmetricTensorTrace (
        MultiNTensorField3D<PRECOMP_T>& tensorField,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedComputeSymmetricTensorTrace (
        MultiNTensorField3D<PRECOMP_T>& tensorField, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void computeVorticity (
        MultiNTensorField3D<PRECOMP_T>& velocity,
        MultiNTensorField3D<PRECOMP_T>& vorticity,
        Box3D domain );

template
void maskedComputeVorticity (
        MultiNTensorField3D<PRECOMP_T>& velocity,
        MultiNTensorField3D<PRECOMP_T>& vorticity, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* computeVorticity (
        MultiNTensorField3D<PRECOMP_T>& velocity,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedComputeVorticity (
        MultiNTensorField3D<PRECOMP_T>& velocity, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void computeBulkVorticity (
        MultiNTensorField3D<PRECOMP_T>& velocity,
        MultiNTensorField3D<PRECOMP_T>& vorticity,
        Box3D domain );

template
void maskedComputeBulkVorticity (
        MultiNTensorField3D<PRECOMP_T>& velocity,
        MultiNTensorField3D<PRECOMP_T>& vorticity, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* computeBulkVorticity (
        MultiNTensorField3D<PRECOMP_T>& velocity,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedComputeBulkVorticity (
        MultiNTensorField3D<PRECOMP_T>& velocity, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void computeStrainRate (
        MultiNTensorField3D<PRECOMP_T>& velocity,
        MultiNTensorField3D<PRECOMP_T>& S,
        Box3D domain );

template
void maskedComputeStrainRate (
        MultiNTensorField3D<PRECOMP_T>& velocity,
        MultiNTensorField3D<PRECOMP_T>& S, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* computeStrainRate (
        MultiNTensorField3D<PRECOMP_T>& velocity,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedComputeStrainRate (
        MultiNTensorField3D<PRECOMP_T>& velocity, MultiNTensorField3D<int>& mask,
        Box3D domain );

template void computeBulkStrainRate (
        MultiNTensorField3D<PRECOMP_T>& velocity,
        MultiNTensorField3D<PRECOMP_T>& S,
        Box3D domain );

template
void maskedComputeBulkStrainRate (
        MultiNTensorField3D<PRECOMP_T>& velocity,
        MultiNTensorField3D<PRECOMP_T>& S, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* computeBulkStrainRate (
        MultiNTensorField3D<PRECOMP_T>& velocity,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedComputeBulkStrainRate (
        MultiNTensorField3D<PRECOMP_T>& velocity, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void add (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<PRECOMP_T>& result,
        Box3D domain );

template
void maskedAdd (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<PRECOMP_T>& result, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* add (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedAdd (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void add (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<PRECOMP_T>& result,
        Box3D domain);

template
void maskedAdd (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<PRECOMP_T>& result, MultiNTensorField3D<int>& mask,
        Box3D domain);

template
MultiNTensorField3D<PRECOMP_T>* add (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedAdd (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
void subtract (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<PRECOMP_T>& result,
        Box3D domain );

template
void maskedSubtract (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<PRECOMP_T>& result, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
MultiNTensorField3D<PRECOMP_T>* subtract (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedSubtract (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
void subtract (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<PRECOMP_T>& result,
        Box3D domain );

template
void maskedSubtract (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<PRECOMP_T>& result, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* subtract (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedSubtract (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void subtract (
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& result,
        Box3D domain );

template
void maskedSubtract (
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& result, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* subtract (
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<PRECOMP_T>& A,
        Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedSubtract (
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<PRECOMP_T>& A, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
void multiply (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<PRECOMP_T>& result,
        Box3D domain );

template
void maskedMultiply (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<PRECOMP_T>& result, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
MultiNTensorField3D<PRECOMP_T>* multiply(MultiNTensorField3D<PRECOMP_T>& A,
                                         MultiNTensorField3D<PRECOMP_T>& B, Box3D domain);

template
MultiNTensorField3D<PRECOMP_T>* maskedMultiply(MultiNTensorField3D<PRECOMP_T>& A,
                                               MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask, Box3D domain);


template
void multiply(MultiNTensorField3D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
              MultiNTensorField3D<PRECOMP_T>& result, Box3D domain);

template
void maskedMultiply(MultiNTensorField3D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
              MultiNTensorField3D<PRECOMP_T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template
MultiNTensorField3D<PRECOMP_T>* multiply(MultiNTensorField3D<PRECOMP_T>& A, PRECOMP_T* alpha, int size, Box3D domain);

template
MultiNTensorField3D<PRECOMP_T>* maskedMultiply(MultiNTensorField3D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
                                               MultiNTensorField3D<int>& mask, Box3D domain);

template
void divide (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<PRECOMP_T>& result,
        Box3D domain );

template
void maskedDivide (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<PRECOMP_T>& result, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
MultiNTensorField3D<PRECOMP_T>* divide(MultiNTensorField3D<PRECOMP_T>& A,
                                       MultiNTensorField3D<PRECOMP_T>& B, Box3D domain);

template
MultiNTensorField3D<PRECOMP_T>* maskedDivide(MultiNTensorField3D<PRECOMP_T>& A,
                                             MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask, Box3D domain);


template
void divide(MultiNTensorField3D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
            MultiNTensorField3D<PRECOMP_T>& result, Box3D domain);

template
void maskedDivide(MultiNTensorField3D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
            MultiNTensorField3D<PRECOMP_T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template
MultiNTensorField3D<PRECOMP_T>* divide(MultiNTensorField3D<PRECOMP_T>& A, PRECOMP_T* alpha, int size, Box3D domain);

template
MultiNTensorField3D<PRECOMP_T>* maskedDivide(MultiNTensorField3D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
                                       MultiNTensorField3D<int>& mask, Box3D domain);


template
void divide(PRECOMP_T* alpha, int size, MultiNTensorField3D<PRECOMP_T>& A,
            MultiNTensorField3D<PRECOMP_T>& result, Box3D domain);

template
void maskedDivide(PRECOMP_T* alpha, int size, MultiNTensorField3D<PRECOMP_T>& A,
            MultiNTensorField3D<PRECOMP_T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template
MultiNTensorField3D<PRECOMP_T>* divide (
      PRECOMP_T* alpha, int size, MultiNTensorField3D<PRECOMP_T>& A, Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedDivide (
      PRECOMP_T* alpha, int size, MultiNTensorField3D<PRECOMP_T>& A, MultiNTensorField3D<int>& mask, Box3D domain );



template
void toThePower (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<PRECOMP_T>& result,
        Box3D domain );

template
void maskedToThePower (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<PRECOMP_T>& result, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
MultiNTensorField3D<PRECOMP_T>* toThePower (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedToThePower (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask, Box3D domain );


template
void toThePower(MultiNTensorField3D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
                MultiNTensorField3D<PRECOMP_T>& result, Box3D domain);

template
void maskedToThePower(MultiNTensorField3D<PRECOMP_T>& A, PRECOMP_T* alpha, int size,
                MultiNTensorField3D<PRECOMP_T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template
MultiNTensorField3D<PRECOMP_T>* toThePower (
        MultiNTensorField3D<PRECOMP_T>& A, PRECOMP_T* alpha, int size, Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedToThePower (
        MultiNTensorField3D<PRECOMP_T>& A, PRECOMP_T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain );


template
void toThePower(PRECOMP_T* alpha, int size, MultiNTensorField3D<PRECOMP_T>& A,
                MultiNTensorField3D<PRECOMP_T>& result, Box3D domain);

template
void maskedToThePower(PRECOMP_T* alpha, int size, MultiNTensorField3D<PRECOMP_T>& A,
                MultiNTensorField3D<PRECOMP_T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template
MultiNTensorField3D<PRECOMP_T>* toThePower (
      PRECOMP_T* alpha, int size, MultiNTensorField3D<PRECOMP_T>& A, Box3D domain );

template
MultiNTensorField3D<PRECOMP_T>* maskedToThePower (
      PRECOMP_T* alpha, int size, MultiNTensorField3D<PRECOMP_T>& A, MultiNTensorField3D<int>& mask, Box3D domain );


template
void equals (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<int>& result,
        Box3D domain );

template
void maskedEquals (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<int>* equals (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
MultiNTensorField3D<int>* maskedEquals (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void equals (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<int>& result,
        Box3D domain);

template
void maskedEquals (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask,
        Box3D domain);

template
MultiNTensorField3D<int>* equals (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box3D domain );

template
MultiNTensorField3D<int>* maskedEquals (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
void lessThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<int>& result,
        Box3D domain );

template
void maskedLessThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<int>* lessThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
MultiNTensorField3D<int>* maskedLessThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void lessThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<int>& result,
        Box3D domain);

template
void maskedLessThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask,
        Box3D domain);

template
MultiNTensorField3D<int>* lessThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box3D domain );

template
MultiNTensorField3D<int>* maskedLessThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
void lessEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<int>& result,
        Box3D domain );

template
void maskedLessEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<int>* lessEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
MultiNTensorField3D<int>* maskedLessEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void lessEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<int>& result,
        Box3D domain);

template
void maskedLessEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask,
        Box3D domain);

template
MultiNTensorField3D<int>* lessEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box3D domain );

template
MultiNTensorField3D<int>* maskedLessEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
void greaterThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<int>& result,
        Box3D domain );

template
void maskedGreaterThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<int>* greaterThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
MultiNTensorField3D<int>* maskedGreaterThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void greaterThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<int>& result,
        Box3D domain);

template
void maskedGreaterThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask,
        Box3D domain);

template
MultiNTensorField3D<int>* greaterThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box3D domain );

template
MultiNTensorField3D<int>* maskedGreaterThan (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
void greaterEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<int>& result,
        Box3D domain );

template
void maskedGreaterEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<int>* greaterEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
MultiNTensorField3D<int>* maskedGreaterEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void greaterEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<int>& result,
        Box3D domain);

template
void maskedGreaterEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask,
        Box3D domain);

template
MultiNTensorField3D<int>* greaterEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box3D domain );

template
MultiNTensorField3D<int>* maskedGreaterEqual (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
void logicalAnd (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<int>& result,
        Box3D domain );

template
void maskedLogicalAnd (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<int>* logicalAnd (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
MultiNTensorField3D<int>* maskedLogicalAnd (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void logicalOr (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<int>& result,
        Box3D domain );

template
void maskedLogicalOr (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
MultiNTensorField3D<int>* logicalOr (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
MultiNTensorField3D<int>* maskedLogicalOr (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void negate (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<int>& result,
        Box3D domain);

template
void maskedNegate (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask,
        Box3D domain);

template
MultiNTensorField3D<int>* negate (
        MultiNTensorField3D<PRECOMP_T>& A,
        Box3D domain );

template
MultiNTensorField3D<int>* maskedNegate (
        MultiNTensorField3D<PRECOMP_T>& A, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
void assignInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
void maskedAssignInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void assignInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box3D domain );

template
void maskedAssignInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void addInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
void maskedAddInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void addInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box3D domain );

template
void maskedAddInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void subtractInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
void maskedSubtractInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void subtractInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box3D domain );

template
void maskedSubtractInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void multiplyInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
void maskedMultiplyInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void multiplyInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box3D domain );

template
void maskedMultiplyInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void divideInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
void maskedDivideInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void divideInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box3D domain );

template
void maskedDivideInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
void toThePowerInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B,
        Box3D domain );

template
void maskedToThePowerInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        MultiNTensorField3D<PRECOMP_T>& B, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void toThePowerInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size,
        Box3D domain );

template
void maskedToThePowerInPlace (
        MultiNTensorField3D<PRECOMP_T>& A,
        PRECOMP_T* alpha, int size, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
PRECOMP_T compute_UPO_ScalarProduct(MultiNTensorField3D<PRECOMP_T>& a, MultiNTensorField3D<PRECOMP_T>& b,
                                    Box3D domain);

template
PRECOMP_T masked_compute_UPO_ScalarProduct(MultiNTensorField3D<PRECOMP_T>& a, MultiNTensorField3D<PRECOMP_T>& b,
                                           MultiNTensorField3D<int>& mask, Box3D domain);

}  // namespace plb
