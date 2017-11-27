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
 * Generator functions for all types of multi-blocks, to make them accessible to SWIG;
 * header file.
 */
#ifndef SWIG_DATA_ANALYSIS_WRAPPER_3D_H
#define SWIG_DATA_ANALYSIS_WRAPPER_3D_H

#include "multiBlock/multiDataField3D.h"

namespace plb {


/* *************** Reductive functions ******************************* */

template<typename T>
void computeAverage(MultiNTensorField3D<T>& vectorField, Box3D domain, T* result, int size);

template<typename T>
void maskedComputeAverage(MultiNTensorField3D<T>& vectorField, MultiNTensorField3D<int>& mask,
                          Box3D domain, T* result, int size);


template<typename T>
void computeMin(MultiNTensorField3D<T>& vectorField, Box3D domain, T* result, int size);

template<typename T>
void maskedComputeMin(MultiNTensorField3D<T>& vectorField, MultiNTensorField3D<int>& mask,
                      Box3D domain, T* result, int size);


template<typename T>
void computeMax(MultiNTensorField3D<T>& vectorField, Box3D domain, T* result, int size);

template<typename T>
void maskedComputeMax(MultiNTensorField3D<T>& vectorField, MultiNTensorField3D<int>& mask,
                      Box3D domain, T* result, int size);


template<typename T>
void computeBoundedAverage(MultiNTensorField3D<T>& vectorField, Box3D domain, T* result, int size);

template<typename T>
void maskedComputeBoundedAverage(MultiNTensorField3D<T>& vectorField, MultiNTensorField3D<int>& mask,
                                 Box3D domain, T* result, int size);


/* *************** Copy-Convert *************************** */

template<typename T1, typename T2>
void maskedCopy( MultiNTensorField3D<T1>& field,
                 MultiNTensorField3D<T2>& convertedField,
                 MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T1, typename T2>
MultiNTensorField3D<T2>* maskedCopyConvert( MultiNTensorField3D<T1>& field, MultiNTensorField3D<int>& mask,
                                            Box3D domain );


/* *************** Component out of a tensor-field ****** */

template<typename T>
void extractComponent(MultiNTensorField3D<T>& tensorField, MultiNTensorField3D<T>& component,
                      Box3D domain, int iComponent);

template<typename T>
void maskedExtractComponent(MultiNTensorField3D<T>& tensorField, MultiNTensorField3D<T>& component,
                            MultiNTensorField3D<int>& mask,
                            Box3D domain, int iComponent);

template<typename T>
MultiNTensorField3D<T>* extractComponent( MultiNTensorField3D<T>& tensorField,
                                          Box3D domain, int iComponent );

template<typename T>
MultiNTensorField3D<T>* maskedExtractComponent( MultiNTensorField3D<T>& tensorField, MultiNTensorField3D<int>& mask,
                                                Box3D domain, int iComponent );


/* *************** Vector-norm of each cell in the field *************** */

template<typename T>
void computeNorm(MultiNTensorField3D<T>& tensorField,
                 MultiNTensorField3D<T>& norm, Box3D domain);

template<typename T>
void maskedComputeNorm(MultiNTensorField3D<T>& tensorField,
                       MultiNTensorField3D<T>& norm,
                       MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* computeNorm(MultiNTensorField3D<T>& tensorField, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedComputeNorm(MultiNTensorField3D<T>& tensorField, MultiNTensorField3D<int>& mask,
                                          Box3D domain);


/* *************** Squared vector-norm of each cell in the field ******** */

template<typename T>
void computeNormSqr(MultiNTensorField3D<T>& tensorField, MultiNTensorField3D<T>& normSqr, Box3D domain);

template<typename T>
void maskedComputeNormSqr(MultiNTensorField3D<T>& tensorField,
                          MultiNTensorField3D<T>& normSqr, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* computeNormSqr(MultiNTensorField3D<T>& tensorField, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedComputeNormSqr(MultiNTensorField3D<T>& tensorField, MultiNTensorField3D<int>& mask,
                                             Box3D domain);


/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template<typename T>
void computeSymmetricTensorNorm(MultiNTensorField3D<T>& tensorField,
                                MultiNTensorField3D<T>& norm,
                                Box3D domain);

template<typename T>
void maskedComputeSymmetricTensorNorm(MultiNTensorField3D<T>& tensorField,
                                      MultiNTensorField3D<T>& norm,
                                      MultiNTensorField3D<int>& mask,
                                      Box3D domain);

template<typename T>
MultiNTensorField3D<T>* computeSymmetricTensorNorm(MultiNTensorField3D<T>& tensorField,
                                                   Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedComputeSymmetricTensorNorm(MultiNTensorField3D<T>& tensorField,
                                                         MultiNTensorField3D<int>& mask,
                                                         Box3D domain);


/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template<typename T>
void computeSymmetricTensorNormSqr(MultiNTensorField3D<T>& tensorField,
                                   MultiNTensorField3D<T>& normSqr,
                                   Box3D domain);

template<typename T>
void maskedComputeSymmetricTensorNormSqr(MultiNTensorField3D<T>& tensorField,
                                         MultiNTensorField3D<T>& normSqr,
                                         MultiNTensorField3D<int>& mask,
                                         Box3D domain);

template<typename T>
MultiNTensorField3D<T>* computeSymmetricTensorNormSqr(MultiNTensorField3D<T>& tensorField,
                                                      Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedComputeSymmetricTensorNormSqr(MultiNTensorField3D<T>& tensorField,
                                                            MultiNTensorField3D<int>& mask,
                                                            Box3D domain);


/* *************** Trace of each symmetric tensor of a field ************ */

template<typename T>
void computeSymmetricTensorTrace(MultiNTensorField3D<T>& tensorField,
                                 MultiNTensorField3D<T>& trace,
                                 Box3D domain);

template<typename T>
void maskedComputeSymmetricTensorTrace(MultiNTensorField3D<T>& tensorField,
                                       MultiNTensorField3D<T>& trace,
                                       MultiNTensorField3D<int>& mask,
                                       Box3D domain);

template<typename T>
MultiNTensorField3D<T>* computeSymmetricTensorTrace(MultiNTensorField3D<T>& tensorField,
                                                    Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedComputeSymmetricTensorTrace(MultiNTensorField3D<T>& tensorField,
                                                          MultiNTensorField3D<int>& mask,
                                                          Box3D domain);


/* *************** Vorticity from Velocity field *********************** */

template<typename T>
void computeVorticity(MultiNTensorField3D<T>& velocity, MultiNTensorField3D<T>& vorticity,
                      Box3D domain);

template<typename T>
void maskedComputeVorticity(MultiNTensorField3D<T>& velocity, MultiNTensorField3D<T>& vorticity,
                            MultiNTensorField3D<int>& mask,
                            Box3D domain);

template<typename T>
MultiNTensorField3D<T>* computeVorticity(MultiNTensorField3D<T>& velocity,
                                         Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedComputeVorticity(MultiNTensorField3D<T>& velocity,
                                              MultiNTensorField3D<int>& mask,
                                              Box3D domain);


/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkVorticity(MultiNTensorField3D<T>& velocity, MultiNTensorField3D<T>& vorticity,
                          Box3D domain);

template<typename T>
void maskedComputeBulkVorticity(MultiNTensorField3D<T>& velocity, MultiNTensorField3D<T>& vorticity,
                                MultiNTensorField3D<int>& mask,
                                Box3D domain);

template<typename T>
MultiNTensorField3D<T>* computeBulkVorticity(MultiNTensorField3D<T>& velocity,
                                             Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedComputeBulkVorticity(MultiNTensorField3D<T>& velocity,
                                                   MultiNTensorField3D<int>& mask,
                                                   Box3D domain);


/* *************** Strain rate from Velocity field ********************* */

template<typename T>
void computeStrainRate(MultiNTensorField3D<T>& velocity,
                       MultiNTensorField3D<T>& S, Box3D domain);

template<typename T>
void maskedComputeStrainRate(MultiNTensorField3D<T>& velocity,
                             MultiNTensorField3D<T>& S,
                             MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* computeStrainRate(MultiNTensorField3D<T>& velocity, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedComputeStrainRate(MultiNTensorField3D<T>& velocity,
                                                MultiNTensorField3D<int>& mask, Box3D domain);


/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkStrainRate(MultiNTensorField3D<T>& velocity, MultiNTensorField3D<T>& S, Box3D domain);

template<typename T>
void maskedComputeBulkStrainRate(MultiNTensorField3D<T>& velocity, MultiNTensorField3D<T>& S,
                                 MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* computeBulkStrainRate(MultiNTensorField3D<T>& velocity, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedComputeBulkStrainRate(MultiNTensorField3D<T>& velocity,
                                                    MultiNTensorField3D<int>& mask, Box3D domain);


/* *************** MultiNTensorField - MultiNTensorField operations *************** */

template<typename T>
void add(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
         MultiNTensorField3D<T>& result, Box3D domain);

template<typename T>
void maskedAdd(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
               MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* add(MultiNTensorField3D<T>& A,
                            MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedAdd(MultiNTensorField3D<T>& A,
                                  MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void add(MultiNTensorField3D<T>& A, T* alpha, int size,
         MultiNTensorField3D<T>& result, Box3D domain);

template<typename T>
void maskedAdd(MultiNTensorField3D<T>& A, T* alpha, int size,
               MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* add(MultiNTensorField3D<T>& A,
                            T* alpha, int size, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedAdd(MultiNTensorField3D<T>& A,
                                  T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void subtract(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
              MultiNTensorField3D<T>& result, Box3D domain);

template<typename T>
void maskedSubtract(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                    MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T> MultiNTensorField3D<T>* subtract(MultiNTensorField3D<T>& A,
                                                      MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedSubtract(MultiNTensorField3D<T>& A,
                                       MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void subtract(MultiNTensorField3D<T>& A, T* alpha, int size,
              MultiNTensorField3D<T>& result, Box3D domain);

template<typename T>
void maskedSubtract(MultiNTensorField3D<T>& A, T* alpha, int size,
                    MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* subtract(MultiNTensorField3D<T>& A,
                                 T* alpha, int size, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedSubtract(MultiNTensorField3D<T>& A,
                                       T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void subtract(T* alpha, int size, MultiNTensorField3D<T>& A,
              MultiNTensorField3D<T>& result, Box3D domain);

template<typename T>
void maskedSubtract(T* alpha, int size, MultiNTensorField3D<T>& A,
                    MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* subtract(T* alpha, int size, MultiNTensorField3D<T>& A, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedSubtract(T* alpha, int size, MultiNTensorField3D<T>& A,
                                       MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void multiply(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
              MultiNTensorField3D<T>& result, Box3D domain);

template<typename T>
void maskedMultiply(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                    MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* multiply(MultiNTensorField3D<T>& A,
                                 MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedMultiply(MultiNTensorField3D<T>& A,
                                       MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void multiply(MultiNTensorField3D<T>& A, T* alpha, int size,
              MultiNTensorField3D<T>& result, Box3D domain);

template<typename T>
void maskedMultiply(MultiNTensorField3D<T>& A, T* alpha, int size,
                    MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* multiply(MultiNTensorField3D<T>& A,
                                 T* alpha, int size, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedMultiply(MultiNTensorField3D<T>& A,
                                       T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void divide(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
            MultiNTensorField3D<T>& result, Box3D domain);

template<typename T>
void maskedDivide(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                  MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* divide(MultiNTensorField3D<T>& A,
                               MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedDivide(MultiNTensorField3D<T>& A,
                                     MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void divide(MultiNTensorField3D<T>& A, T* alpha, int size,
            MultiNTensorField3D<T>& result, Box3D domain);

template<typename T>
void maskedDivide(MultiNTensorField3D<T>& A, T* alpha, int size,
                  MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* divide(MultiNTensorField3D<T>& A,
                               T* alpha, int size, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedDivide(MultiNTensorField3D<T>& A,
                                     T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void divide(T* alpha, int size, MultiNTensorField3D<T>& A,
            MultiNTensorField3D<T>& result, Box3D domain);

template<typename T>
void maskedDivide(T* alpha, int size, MultiNTensorField3D<T>& A,
                  MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* divide(T* alpha, int size,
                               MultiNTensorField3D<T>& A, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedDivide(T* alpha, int size,
                                     MultiNTensorField3D<T>& A, MultiNTensorField3D<int>& mask, Box3D domain);



template<typename T>
void toThePower(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                MultiNTensorField3D<T>& result, Box3D domain);

template<typename T>
void maskedToThePower(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                      MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* toThePower(MultiNTensorField3D<T>& A,
                                   MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedToThePower(MultiNTensorField3D<T>& A,
                                         MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void toThePower(MultiNTensorField3D<T>& A, T* alpha, int size,
                MultiNTensorField3D<T>& result, Box3D domain);

template<typename T>
void maskedToThePower(MultiNTensorField3D<T>& A, T* alpha, int size,
                      MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* toThePower(MultiNTensorField3D<T>& A,
                                   T* alpha, int size, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedToThePower(MultiNTensorField3D<T>& A,
                                         T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void toThePower(T* alpha, int size, MultiNTensorField3D<T>& A,
                MultiNTensorField3D<T>& result, Box3D domain);

template<typename T>
void maskedToThePower(T* alpha, int size, MultiNTensorField3D<T>& A,
                      MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* toThePower(T* alpha, int size,
                                   MultiNTensorField3D<T>& A, Box3D domain);

template<typename T>
MultiNTensorField3D<T>* maskedToThePower(T* alpha, int size,
                                         MultiNTensorField3D<T>& A, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void equals(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
            MultiNTensorField3D<int>& result, Box3D domain);

template<typename T>
void maskedEquals(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                  MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* equals(MultiNTensorField3D<T>& A,
                                 MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* maskedEquals(MultiNTensorField3D<T>& A,
                                       MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
void equals(MultiNTensorField3D<T>& A, T* alpha, int size,
            MultiNTensorField3D<int>& result, Box3D domain);

template<typename T>
void maskedEquals(MultiNTensorField3D<T>& A, T* alpha, int size,
                  MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* equals(MultiNTensorField3D<T>& A,
                                 T* alpha, int size, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* maskedEquals(MultiNTensorField3D<T>& A,
                                       T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void lessThan(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
              MultiNTensorField3D<int>& result, Box3D domain);

template<typename T>
void maskedLessThan(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                    MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* lessThan(MultiNTensorField3D<T>& A,
                                   MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* maskedLessThan(MultiNTensorField3D<T>& A,
                                         MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
void lessThan(MultiNTensorField3D<T>& A, T* alpha, int size,
              MultiNTensorField3D<int>& result, Box3D domain);

template<typename T>
void maskedLessThan(MultiNTensorField3D<T>& A, T* alpha, int size,
                    MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* lessThan(MultiNTensorField3D<T>& A,
                                   T* alpha, int size, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* maskedLessThan(MultiNTensorField3D<T>& A,
                                         T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void lessEqual(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
               MultiNTensorField3D<int>& result, Box3D domain);

template<typename T>
void maskedLessEqual(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                     MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* lessEqual(MultiNTensorField3D<T>& A,
                                    MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* maskedLessEqual(MultiNTensorField3D<T>& A,
                                          MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
void lessEqual(MultiNTensorField3D<T>& A, T* alpha, int size,
               MultiNTensorField3D<int>& result, Box3D domain);

template<typename T>
void maskedLessEqual(MultiNTensorField3D<T>& A, T* alpha, int size,
                     MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* lessEqual(MultiNTensorField3D<T>& A,
                                    T* alpha, int size, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* maskedLessEqual(MultiNTensorField3D<T>& A,
                                          T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void greaterThan(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                 MultiNTensorField3D<int>& result, Box3D domain);

template<typename T>
void maskedGreaterThan(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                       MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* greaterThan(MultiNTensorField3D<T>& A,
                                      MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* maskedGreaterThan(MultiNTensorField3D<T>& A,
                                            MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
void greaterThan(MultiNTensorField3D<T>& A, T* alpha, int size,
                 MultiNTensorField3D<int>& result, Box3D domain);

template<typename T>
void maskedGreaterThan(MultiNTensorField3D<T>& A, T* alpha, int size,
                       MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* greaterThan(MultiNTensorField3D<T>& A,
                                      T* alpha, int size, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* maskedGreaterThan(MultiNTensorField3D<T>& A,
                                            T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);



template<typename T>
void greaterEqual(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                  MultiNTensorField3D<int>& result, Box3D domain);

template<typename T>
void maskedGreaterEqual(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* greaterEqual(MultiNTensorField3D<T>& A,
                                       MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* maskedGreaterEqual(MultiNTensorField3D<T>& A,
                                             MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
void greaterEqual(MultiNTensorField3D<T>& A, T* alpha, int size,
                  MultiNTensorField3D<int>& result, Box3D domain);

template<typename T>
void maskedGreaterEqual(MultiNTensorField3D<T>& A, T* alpha, int size,
                        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* greaterEqual(MultiNTensorField3D<T>& A,
                                       T* alpha, int size, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* maskedGreaterEqual(MultiNTensorField3D<T>& A,
                                             T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void logicalAnd(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                MultiNTensorField3D<int>& result, Box3D domain);

template<typename T>
void maskedLogicalAnd(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                      MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* logicalAnd(MultiNTensorField3D<T>& A,
                                     MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* maskedLogicalAnd(MultiNTensorField3D<T>& A,
                                           MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void logicalOr(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
               MultiNTensorField3D<int>& result, Box3D domain);

template<typename T>
void maskedLogicalOr(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                     MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* logicalOr(MultiNTensorField3D<T>& A,
                                    MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* maskedLogicalOr(MultiNTensorField3D<T>& A,
                                          MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);


template<typename T>
void negate(MultiNTensorField3D<T>& A,
            MultiNTensorField3D<int>& result, Box3D domain);

template<typename T>
void maskedNegate(MultiNTensorField3D<T>& A,
                  MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T>
MultiNTensorField3D<int>* negate(MultiNTensorField3D<T>& A,
                                 Box3D domain);

template<typename T>
MultiNTensorField3D<int>* maskedNegate(MultiNTensorField3D<T>& A, MultiNTensorField3D<int>& mask,
                                       Box3D domain);


/* *************** MultiNTensorField - MultiNTensorField inplace operations *************** */

template<typename T>
void assignInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
void maskedAssignInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);
      
template<typename T>
void assignInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, Box3D domain);

template<typename T>
void maskedAssignInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);
      

template<typename T>
void addInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
void maskedAddInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);
      
template<typename T>
void addInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, Box3D domain);

template<typename T>
void maskedAddInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);
      

template<typename T>
void subtractInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
void maskedSubtractInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);
      
template<typename T>
void subtractInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, Box3D domain);

template<typename T>
void maskedSubtractInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);
      

template<typename T>
void multiplyInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
void maskedMultiplyInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);
      
template<typename T>
void multiplyInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, Box3D domain);

template<typename T>
void maskedMultiplyInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);
      

template<typename T>
void divideInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
void maskedDivideInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);
      
template<typename T>
void divideInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, Box3D domain);

template<typename T>
void maskedDivideInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);

      
template<typename T>
void toThePowerInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, Box3D domain);

template<typename T>
void maskedToThePowerInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain);
      
template<typename T>
void toThePowerInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, Box3D domain);

template<typename T>
void maskedToThePowerInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain);
      
/* *************** UPO ******************************* */

template<typename T>
T compute_UPO_ScalarProduct(MultiNTensorField3D<T>& a, MultiNTensorField3D<T>& b,
                            Box3D domain);

template<typename T>
T masked_compute_UPO_ScalarProduct(MultiNTensorField3D<T>& a, MultiNTensorField3D<T>& b,
                                   MultiNTensorField3D<int>& mask, Box3D domain);

}  // namespace plb

#endif  // SWIG_DATA_ANALYSIS_WRAPPER_3D_H
