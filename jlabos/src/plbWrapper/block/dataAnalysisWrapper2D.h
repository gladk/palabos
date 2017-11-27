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
#ifndef SWIG_DATA_ANALYSIS_WRAPPER_2D_H
#define SWIG_DATA_ANALYSIS_WRAPPER_2D_H

#include "multiBlock/multiDataField2D.h"

namespace plb {


/* *************** Reductive functions ******************************* */

template<typename T>
void computeAverage(MultiNTensorField2D<T>& vectorField, Box2D domain, T* result, int size);

template<typename T>
void maskedComputeAverage(MultiNTensorField2D<T>& vectorField, MultiNTensorField2D<int>& mask,
                          Box2D domain, T* result, int size);


template<typename T>
void computeMin(MultiNTensorField2D<T>& vectorField, Box2D domain, T* result, int size);

template<typename T>
void maskedComputeMin(MultiNTensorField2D<T>& vectorField, MultiNTensorField2D<int>& mask,
                      Box2D domain, T* result, int size);


template<typename T>
void computeMax(MultiNTensorField2D<T>& vectorField, Box2D domain, T* result, int size);

template<typename T>
void maskedComputeMax(MultiNTensorField2D<T>& vectorField, MultiNTensorField2D<int>& mask,
                      Box2D domain, T* result, int size);


template<typename T>
void computeBoundedAverage(MultiNTensorField2D<T>& vectorField, Box2D domain, T* result, int size);

template<typename T>
void maskedComputeBoundedAverage(MultiNTensorField2D<T>& vectorField, MultiNTensorField2D<int>& mask,
                                 Box2D domain, T* result, int size);


/* *************** Copy-Convert *************************** */

template<typename T1, typename T2>
void maskedCopy( MultiNTensorField2D<T1>& field,
                 MultiNTensorField2D<T2>& convertedField,
                 MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T1, typename T2>
MultiNTensorField2D<T2>* maskedCopyConvert( MultiNTensorField2D<T1>& field, MultiNTensorField2D<int>& mask,
                                            Box2D domain );


/* *************** Component out of a tensor-field ****** */

template<typename T>
void extractComponent(MultiNTensorField2D<T>& tensorField, MultiNTensorField2D<T>& component,
                      Box2D domain, int iComponent);

template<typename T>
void maskedExtractComponent(MultiNTensorField2D<T>& tensorField, MultiNTensorField2D<T>& component,
                            MultiNTensorField2D<int>& mask,
                            Box2D domain, int iComponent);

template<typename T>
MultiNTensorField2D<T>* extractComponent( MultiNTensorField2D<T>& tensorField,
                                          Box2D domain, int iComponent );

template<typename T>
MultiNTensorField2D<T>* maskedExtractComponent( MultiNTensorField2D<T>& tensorField, MultiNTensorField2D<int>& mask,
                                                Box2D domain, int iComponent );


/* *************** Vector-norm of each cell in the field *************** */

template<typename T>
void computeNorm(MultiNTensorField2D<T>& tensorField,
                 MultiNTensorField2D<T>& norm, Box2D domain);

template<typename T>
void maskedComputeNorm(MultiNTensorField2D<T>& tensorField,
                       MultiNTensorField2D<T>& norm,
                       MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* computeNorm(MultiNTensorField2D<T>& tensorField, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedComputeNorm(MultiNTensorField2D<T>& tensorField, MultiNTensorField2D<int>& mask,
                                          Box2D domain);


/* *************** Squared vector-norm of each cell in the field ******** */

template<typename T>
void computeNormSqr(MultiNTensorField2D<T>& tensorField, MultiNTensorField2D<T>& normSqr, Box2D domain);

template<typename T>
void maskedComputeNormSqr(MultiNTensorField2D<T>& tensorField,
                          MultiNTensorField2D<T>& normSqr, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* computeNormSqr(MultiNTensorField2D<T>& tensorField, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedComputeNormSqr(MultiNTensorField2D<T>& tensorField, MultiNTensorField2D<int>& mask,
                                             Box2D domain);


/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template<typename T>
void computeSymmetricTensorNorm(MultiNTensorField2D<T>& tensorField,
                                MultiNTensorField2D<T>& norm,
                                Box2D domain);

template<typename T>
void maskedComputeSymmetricTensorNorm(MultiNTensorField2D<T>& tensorField,
                                      MultiNTensorField2D<T>& norm,
                                      MultiNTensorField2D<int>& mask,
                                      Box2D domain);

template<typename T>
MultiNTensorField2D<T>* computeSymmetricTensorNorm(MultiNTensorField2D<T>& tensorField,
                                                   Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedComputeSymmetricTensorNorm(MultiNTensorField2D<T>& tensorField,
                                                         MultiNTensorField2D<int>& mask,
                                                         Box2D domain);


/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template<typename T>
void computeSymmetricTensorNormSqr(MultiNTensorField2D<T>& tensorField,
                                   MultiNTensorField2D<T>& normSqr,
                                   Box2D domain);

template<typename T>
void maskedComputeSymmetricTensorNormSqr(MultiNTensorField2D<T>& tensorField,
                                         MultiNTensorField2D<T>& normSqr,
                                         MultiNTensorField2D<int>& mask,
                                         Box2D domain);

template<typename T>
MultiNTensorField2D<T>* computeSymmetricTensorNormSqr(MultiNTensorField2D<T>& tensorField,
                                                      Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedComputeSymmetricTensorNormSqr(MultiNTensorField2D<T>& tensorField,
                                                            MultiNTensorField2D<int>& mask,
                                                            Box2D domain);


/* *************** Trace of each symmetric tensor of a field ************ */

template<typename T>
void computeSymmetricTensorTrace(MultiNTensorField2D<T>& tensorField,
                                 MultiNTensorField2D<T>& trace,
                                 Box2D domain);

template<typename T>
void maskedComputeSymmetricTensorTrace(MultiNTensorField2D<T>& tensorField,
                                       MultiNTensorField2D<T>& trace,
                                       MultiNTensorField2D<int>& mask,
                                       Box2D domain);

template<typename T>
MultiNTensorField2D<T>* computeSymmetricTensorTrace(MultiNTensorField2D<T>& tensorField,
                                                    Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedComputeSymmetricTensorTrace(MultiNTensorField2D<T>& tensorField,
                                                          MultiNTensorField2D<int>& mask,
                                                          Box2D domain);


/* *************** Vorticity from Velocity field *********************** */

template<typename T>
void computeVorticity(MultiNTensorField2D<T>& velocity, MultiNTensorField2D<T>& vorticity,
                      Box2D domain);

template<typename T>
void maskedComputeVorticity(MultiNTensorField2D<T>& velocity, MultiNTensorField2D<T>& vorticity,
                            MultiNTensorField2D<int>& mask,
                            Box2D domain);

template<typename T>
MultiNTensorField2D<T>* computeVorticity(MultiNTensorField2D<T>& velocity,
                                         Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedComputeVorticity(MultiNTensorField2D<T>& velocity,
                                              MultiNTensorField2D<int>& mask,
                                              Box2D domain);


/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkVorticity(MultiNTensorField2D<T>& velocity, MultiNTensorField2D<T>& vorticity,
                          Box2D domain);

template<typename T>
void maskedComputeBulkVorticity(MultiNTensorField2D<T>& velocity, MultiNTensorField2D<T>& vorticity,
                                MultiNTensorField2D<int>& mask,
                                Box2D domain);

template<typename T>
MultiNTensorField2D<T>* computeBulkVorticity(MultiNTensorField2D<T>& velocity,
                                             Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedComputeBulkVorticity(MultiNTensorField2D<T>& velocity,
                                                   MultiNTensorField2D<int>& mask,
                                                   Box2D domain);


/* *************** Strain rate from Velocity field ********************* */

template<typename T>
void computeStrainRate(MultiNTensorField2D<T>& velocity,
                       MultiNTensorField2D<T>& S, Box2D domain);

template<typename T>
void maskedComputeStrainRate(MultiNTensorField2D<T>& velocity,
                             MultiNTensorField2D<T>& S,
                             MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* computeStrainRate(MultiNTensorField2D<T>& velocity, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedComputeStrainRate(MultiNTensorField2D<T>& velocity,
                                                MultiNTensorField2D<int>& mask, Box2D domain);


/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkStrainRate(MultiNTensorField2D<T>& velocity, MultiNTensorField2D<T>& S, Box2D domain);

template<typename T>
void maskedComputeBulkStrainRate(MultiNTensorField2D<T>& velocity, MultiNTensorField2D<T>& S,
                                 MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* computeBulkStrainRate(MultiNTensorField2D<T>& velocity, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedComputeBulkStrainRate(MultiNTensorField2D<T>& velocity,
                                                    MultiNTensorField2D<int>& mask, Box2D domain);


/* *************** MultiNTensorField - MultiNTensorField operations *************** */

template<typename T>
void add(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
         MultiNTensorField2D<T>& result, Box2D domain);

template<typename T>
void maskedAdd(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
               MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* add(MultiNTensorField2D<T>& A,
                            MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedAdd(MultiNTensorField2D<T>& A,
                                  MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void add(MultiNTensorField2D<T>& A, T* alpha, int size,
         MultiNTensorField2D<T>& result, Box2D domain);

template<typename T>
void maskedAdd(MultiNTensorField2D<T>& A, T* alpha, int size,
               MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* add(MultiNTensorField2D<T>& A,
                            T* alpha, int size, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedAdd(MultiNTensorField2D<T>& A,
                                  T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void subtract(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
              MultiNTensorField2D<T>& result, Box2D domain);

template<typename T>
void maskedSubtract(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                    MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T> MultiNTensorField2D<T>* subtract(MultiNTensorField2D<T>& A,
                                                      MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedSubtract(MultiNTensorField2D<T>& A,
                                       MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void subtract(MultiNTensorField2D<T>& A, T* alpha, int size,
              MultiNTensorField2D<T>& result, Box2D domain);

template<typename T>
void maskedSubtract(MultiNTensorField2D<T>& A, T* alpha, int size,
                    MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* subtract(MultiNTensorField2D<T>& A,
                                 T* alpha, int size, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedSubtract(MultiNTensorField2D<T>& A,
                                       T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void subtract(T* alpha, int size, MultiNTensorField2D<T>& A,
              MultiNTensorField2D<T>& result, Box2D domain);

template<typename T>
void maskedSubtract(T* alpha, int size, MultiNTensorField2D<T>& A,
                    MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* subtract(T* alpha, int size, MultiNTensorField2D<T>& A, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedSubtract(T* alpha, int size, MultiNTensorField2D<T>& A,
                                       MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void multiply(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
              MultiNTensorField2D<T>& result, Box2D domain);

template<typename T>
void maskedMultiply(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                    MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* multiply(MultiNTensorField2D<T>& A,
                                 MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedMultiply(MultiNTensorField2D<T>& A,
                                       MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void multiply(MultiNTensorField2D<T>& A, T* alpha, int size,
              MultiNTensorField2D<T>& result, Box2D domain);

template<typename T>
void maskedMultiply(MultiNTensorField2D<T>& A, T* alpha, int size,
                    MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* multiply(MultiNTensorField2D<T>& A,
                                 T* alpha, int size, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedMultiply(MultiNTensorField2D<T>& A,
                                       T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void divide(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
            MultiNTensorField2D<T>& result, Box2D domain);

template<typename T>
void maskedDivide(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                  MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* divide(MultiNTensorField2D<T>& A,
                               MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedDivide(MultiNTensorField2D<T>& A,
                                     MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void divide(MultiNTensorField2D<T>& A, T* alpha, int size,
            MultiNTensorField2D<T>& result, Box2D domain);

template<typename T>
void maskedDivide(MultiNTensorField2D<T>& A, T* alpha, int size,
                  MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* divide(MultiNTensorField2D<T>& A,
                               T* alpha, int size, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedDivide(MultiNTensorField2D<T>& A,
                                     T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void divide(T* alpha, int size, MultiNTensorField2D<T>& A,
            MultiNTensorField2D<T>& result, Box2D domain);

template<typename T>
void maskedDivide(T* alpha, int size, MultiNTensorField2D<T>& A,
                  MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* divide(T* alpha, int size,
                               MultiNTensorField2D<T>& A, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedDivide(T* alpha, int size,
                                     MultiNTensorField2D<T>& A, MultiNTensorField2D<int>& mask, Box2D domain);



template<typename T>
void toThePower(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                MultiNTensorField2D<T>& result, Box2D domain);

template<typename T>
void maskedToThePower(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                      MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* toThePower(MultiNTensorField2D<T>& A,
                                   MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedToThePower(MultiNTensorField2D<T>& A,
                                         MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void toThePower(MultiNTensorField2D<T>& A, T* alpha, int size,
                MultiNTensorField2D<T>& result, Box2D domain);

template<typename T>
void maskedToThePower(MultiNTensorField2D<T>& A, T* alpha, int size,
                      MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* toThePower(MultiNTensorField2D<T>& A,
                                   T* alpha, int size, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedToThePower(MultiNTensorField2D<T>& A,
                                         T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void toThePower(T* alpha, int size, MultiNTensorField2D<T>& A,
                MultiNTensorField2D<T>& result, Box2D domain);

template<typename T>
void maskedToThePower(T* alpha, int size, MultiNTensorField2D<T>& A,
                      MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* toThePower(T* alpha, int size,
                                   MultiNTensorField2D<T>& A, Box2D domain);

template<typename T>
MultiNTensorField2D<T>* maskedToThePower(T* alpha, int size,
                                         MultiNTensorField2D<T>& A, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void equals(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
            MultiNTensorField2D<int>& result, Box2D domain);

template<typename T>
void maskedEquals(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                  MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* equals(MultiNTensorField2D<T>& A,
                                 MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* maskedEquals(MultiNTensorField2D<T>& A,
                                       MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
void equals(MultiNTensorField2D<T>& A, T* alpha, int size,
            MultiNTensorField2D<int>& result, Box2D domain);

template<typename T>
void maskedEquals(MultiNTensorField2D<T>& A, T* alpha, int size,
                  MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* equals(MultiNTensorField2D<T>& A,
                                 T* alpha, int size, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* maskedEquals(MultiNTensorField2D<T>& A,
                                       T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void lessThan(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
              MultiNTensorField2D<int>& result, Box2D domain);

template<typename T>
void maskedLessThan(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                    MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* lessThan(MultiNTensorField2D<T>& A,
                                   MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* maskedLessThan(MultiNTensorField2D<T>& A,
                                         MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
void lessThan(MultiNTensorField2D<T>& A, T* alpha, int size,
              MultiNTensorField2D<int>& result, Box2D domain);

template<typename T>
void maskedLessThan(MultiNTensorField2D<T>& A, T* alpha, int size,
                    MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* lessThan(MultiNTensorField2D<T>& A,
                                   T* alpha, int size, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* maskedLessThan(MultiNTensorField2D<T>& A,
                                         T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void lessEqual(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
               MultiNTensorField2D<int>& result, Box2D domain);

template<typename T>
void maskedLessEqual(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                     MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* lessEqual(MultiNTensorField2D<T>& A,
                                    MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* maskedLessEqual(MultiNTensorField2D<T>& A,
                                          MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
void lessEqual(MultiNTensorField2D<T>& A, T* alpha, int size,
               MultiNTensorField2D<int>& result, Box2D domain);

template<typename T>
void maskedLessEqual(MultiNTensorField2D<T>& A, T* alpha, int size,
                     MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* lessEqual(MultiNTensorField2D<T>& A,
                                    T* alpha, int size, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* maskedLessEqual(MultiNTensorField2D<T>& A,
                                          T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void greaterThan(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                 MultiNTensorField2D<int>& result, Box2D domain);

template<typename T>
void maskedGreaterThan(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                       MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* greaterThan(MultiNTensorField2D<T>& A,
                                      MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* maskedGreaterThan(MultiNTensorField2D<T>& A,
                                            MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
void greaterThan(MultiNTensorField2D<T>& A, T* alpha, int size,
                 MultiNTensorField2D<int>& result, Box2D domain);

template<typename T>
void maskedGreaterThan(MultiNTensorField2D<T>& A, T* alpha, int size,
                       MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* greaterThan(MultiNTensorField2D<T>& A,
                                      T* alpha, int size, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* maskedGreaterThan(MultiNTensorField2D<T>& A,
                                            T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);



template<typename T>
void greaterEqual(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                  MultiNTensorField2D<int>& result, Box2D domain);

template<typename T>
void maskedGreaterEqual(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* greaterEqual(MultiNTensorField2D<T>& A,
                                       MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* maskedGreaterEqual(MultiNTensorField2D<T>& A,
                                             MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
void greaterEqual(MultiNTensorField2D<T>& A, T* alpha, int size,
                  MultiNTensorField2D<int>& result, Box2D domain);

template<typename T>
void maskedGreaterEqual(MultiNTensorField2D<T>& A, T* alpha, int size,
                        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* greaterEqual(MultiNTensorField2D<T>& A,
                                       T* alpha, int size, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* maskedGreaterEqual(MultiNTensorField2D<T>& A,
                                             T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void logicalAnd(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                MultiNTensorField2D<int>& result, Box2D domain);

template<typename T>
void maskedLogicalAnd(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                      MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* logicalAnd(MultiNTensorField2D<T>& A,
                                     MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* maskedLogicalAnd(MultiNTensorField2D<T>& A,
                                           MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void logicalOr(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
               MultiNTensorField2D<int>& result, Box2D domain);

template<typename T>
void maskedLogicalOr(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                     MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* logicalOr(MultiNTensorField2D<T>& A,
                                    MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* maskedLogicalOr(MultiNTensorField2D<T>& A,
                                          MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);


template<typename T>
void negate(MultiNTensorField2D<T>& A,
            MultiNTensorField2D<int>& result, Box2D domain);

template<typename T>
void maskedNegate(MultiNTensorField2D<T>& A,
                  MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T>
MultiNTensorField2D<int>* negate(MultiNTensorField2D<T>& A,
                                 Box2D domain);

template<typename T>
MultiNTensorField2D<int>* maskedNegate(MultiNTensorField2D<T>& A, MultiNTensorField2D<int>& mask,
                                       Box2D domain);


/* *************** MultiNTensorField - MultiNTensorField inplace operations *************** */

template<typename T>
void assignInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
void maskedAssignInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);
      
template<typename T>
void assignInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, Box2D domain);

template<typename T>
void maskedAssignInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);
      

template<typename T>
void addInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
void maskedAddInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);
      
template<typename T>
void addInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, Box2D domain);

template<typename T>
void maskedAddInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);
      

template<typename T>
void subtractInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
void maskedSubtractInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);
      
template<typename T>
void subtractInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, Box2D domain);

template<typename T>
void maskedSubtractInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);
      

template<typename T>
void multiplyInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
void maskedMultiplyInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);
      
template<typename T>
void multiplyInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, Box2D domain);

template<typename T>
void maskedMultiplyInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);
      

template<typename T>
void divideInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
void maskedDivideInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);
      
template<typename T>
void divideInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, Box2D domain);

template<typename T>
void maskedDivideInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);

      
template<typename T>
void toThePowerInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, Box2D domain);

template<typename T>
void maskedToThePowerInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain);
      
template<typename T>
void toThePowerInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, Box2D domain);

template<typename T>
void maskedToThePowerInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain);

      
}  // namespace plb

#endif  // SWIG_DATA_ANALYSIS_WRAPPER_2D_H
