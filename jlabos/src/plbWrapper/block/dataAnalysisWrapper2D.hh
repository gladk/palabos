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

#ifndef SWIG_DATA_ANALYSIS_WRAPPER_2D_HH
#define SWIG_DATA_ANALYSIS_WRAPPER_2D_HH

#include "plbWrapper/block/dataAnalysisWrapper2D.h"
#include "plbWrapper/block/dataAnalysisFunctional2D.h"
#include "plbWrapper/block/dataInitializerFunctional2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "multiBlock/multiBlockGenerator2D.h"

namespace plb {


/* *************** Reductive functions ******************************* */

template<typename T>
void computeAverage(MultiNTensorField2D<T>& vectorField, Box2D domain, T* result, int size)
{
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    BoxNTensorSumFunctional2D<T> functional(vectorField.getNdim());
    applyProcessingFunctional(functional, domain, vectorField);
    std::vector<T> const& sumVector = functional.getSumVector();
    for (plint iDim=0; iDim<vectorField.getNdim(); ++iDim) {
        result[iDim] = sumVector[iDim] / (T) domain.nCells();
    }
}

template<typename T>
void maskedComputeAverage(MultiNTensorField2D<T>& vectorField, MultiNTensorField2D<int>& mask, Box2D domain, T* result, int size)
{
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    MaskedBoxNTensorSumFunctional2D<T> functional(vectorField.getNdim());
    applyProcessingFunctional(functional, domain, vectorField, mask);
    std::vector<T> const& sumVector = functional.getSumVector();
    for (plint iDim=0; iDim<vectorField.getNdim(); ++iDim) {
        result[iDim] = sumVector[iDim] / (T) domain.nCells();
    }
}

template<typename T>
void computeMin(MultiNTensorField2D<T>& vectorField, Box2D domain, T* result, int size)
{
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    BoxNTensorMinFunctional2D<T> functional(vectorField.getNdim());
    applyProcessingFunctional(functional, domain, vectorField);
    std::vector<T> const& minVector = functional.getMinVector();
    for (plint iDim=0; iDim<vectorField.getNdim(); ++iDim) {
        result[iDim] = minVector[iDim];
    }
}

template<typename T>
void maskedComputeMin(MultiNTensorField2D<T>& vectorField, MultiNTensorField2D<int>& mask, Box2D domain, T* result, int size)
{
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    MaskedBoxNTensorMinFunctional2D<T> functional(vectorField.getNdim());
    applyProcessingFunctional(functional, domain, vectorField, mask);
    std::vector<T> const& minVector = functional.getMinVector();
    for (plint iDim=0; iDim<vectorField.getNdim(); ++iDim) {
        result[iDim] = minVector[iDim];
    }
}


template<typename T>
void computeMax(MultiNTensorField2D<T>& vectorField, Box2D domain, T* result, int size)
{
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    BoxNTensorMaxFunctional2D<T> functional(vectorField.getNdim());
    applyProcessingFunctional(functional, domain, vectorField);
    std::vector<T> const& maxVector = functional.getMaxVector();
    for (plint iDim=0; iDim<vectorField.getNdim(); ++iDim) {
        result[iDim] = maxVector[iDim];
    }
}

template<typename T>
void maskedComputeMax(MultiNTensorField2D<T>& vectorField, MultiNTensorField2D<int>& mask, Box2D domain, T* result, int size)
{
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    MaskedBoxNTensorMaxFunctional2D<T> functional(vectorField.getNdim());
    applyProcessingFunctional(functional, domain, vectorField, mask);
    std::vector<T> const& maxVector = functional.getMaxVector();
    for (plint iDim=0; iDim<vectorField.getNdim(); ++iDim) {
        result[iDim] = maxVector[iDim];
    }
}

template<typename T>
void computeBoundedAverage(MultiNTensorField2D<T>& vectorField, Box2D domain, T* result, int size) {
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    BoundedBoxNTensorSumFunctional2D<T> functional(vectorField.getNdim());
    plint envelopeWidth=1;
    applyProcessingFunctional(functional, domain, vectorField, envelopeWidth);
    std::vector<T> sumVector(functional.getSumVector());
    for (plint iDim=0; iDim<vectorField.getNdim(); ++iDim) {
        result[iDim] = sumVector[iDim] / 
             (T) ( (domain.getNx()-1)*(domain.getNy()-1) );
    }
}

template<typename T>
void maskedComputeBoundedAverage (
        MultiNTensorField2D<T>& vectorField, MultiNTensorField2D<int>& mask,
        Box2D domain, T* result, int size)
{
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    BoundedMaskedBoxNTensorSumFunctional2D<T> functional(vectorField.getNdim());
    plint envelopeWidth=1;
    applyProcessingFunctional(functional, domain, vectorField, mask, envelopeWidth);
    std::vector<T> sumVector(functional.getSumVector());
    for (plint iDim=0; iDim<vectorField.getNdim(); ++iDim) {
        result[iDim] = sumVector[iDim] / 
             (T) ( (domain.getNx()-1)*(domain.getNy()-1) );
    }
}


/* *************** Copy-Convert *************************** */

template<typename T1, typename T2>
void maskedCopy( MultiNTensorField2D<T1>& field,
                 MultiNTensorField2D<T2>& convertedField, MultiNTensorField2D<int>& mask,
                 Box2D domain)
{
    applyProcessingFunctional (
            new MaskedCopyConvertNTensorFunctional2D<T1,T2>, domain, field, convertedField, mask );
}


template<typename T1, typename T2>
MultiNTensorField2D<T2>* maskedCopyConvert( MultiNTensorField2D<T1>& field, MultiNTensorField2D<int>& mask,
                                            Box2D domain)
{
    MultiNTensorField2D<T2>* convertedField
        = generateMultiNTensorField<T2>(field, domain, field.getNdim());
    maskedCopy(field, *convertedField, mask, domain);
    return convertedField;
}

/* *************** Component (scalar-field) out of a tensor-field ****** */

template<typename T>
void extractComponent(MultiNTensorField2D<T>& vectorField, MultiNTensorField2D<T>& component,
                      Box2D domain, int iComponent)
{
    applyProcessingFunctional (
            new ExtractNTensorComponentFunctional2D<T>(iComponent),
            domain, component, vectorField );
}

template<typename T>
void maskedExtractComponent(MultiNTensorField2D<T>& vectorField, MultiNTensorField2D<T>& component, MultiNTensorField2D<int>& mask,
                            Box2D domain, int iComponent)
{
    applyProcessingFunctional (
            new MaskedExtractNTensorComponentFunctional2D<T>(iComponent),
            domain, component, vectorField, mask );
}

template<typename T>
MultiNTensorField2D<T>* extractComponent (
              MultiNTensorField2D<T>& vectorField,
              Box2D domain, int iComponent )
{
    MultiNTensorField2D<T>* component
        = generateMultiNTensorField<T>(vectorField, domain, 1);
    extractComponent(vectorField, *component, domain, iComponent);
    return component;
}

template<typename T>
MultiNTensorField2D<T>* maskedExtractComponent (
        MultiNTensorField2D<T>& vectorField, MultiNTensorField2D<int>& mask,
        Box2D domain, int iComponent )
{
    MultiNTensorField2D<T>* component
        = generateMultiNTensorField<T>(vectorField, domain, 1);
    maskedExtractComponent(vectorField, *component, mask, domain, iComponent);
    return component;
}

/* *************** Vector-norm of each cell in the field *************** */

template<typename T>
void computeNorm(MultiNTensorField2D<T>& vectorField,
                 MultiNTensorField2D<T>& norm, Box2D domain)
{
    applyProcessingFunctional (
            new ComputeNTensorNormFunctional2D<T>, domain, norm, vectorField );
}

template<typename T>
void maskedComputeNorm(MultiNTensorField2D<T>& vectorField,
                 MultiNTensorField2D<T>& norm, MultiNTensorField2D<int>& mask, Box2D domain)
{
    applyProcessingFunctional (
            new MaskedComputeNTensorNormFunctional2D<T>, domain, norm, vectorField, mask );
}

template<typename T>
MultiNTensorField2D<T>* computeNorm(MultiNTensorField2D<T>& vectorField, Box2D domain)
{
    MultiNTensorField2D<T>* norm
        = generateMultiNTensorField<T>(vectorField, domain, 1);
    computeNorm(vectorField, *norm, domain);
    return norm;
}

template<typename T>
MultiNTensorField2D<T>* maskedComputeNorm(MultiNTensorField2D<T>& vectorField, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<T>* norm
        = generateMultiNTensorField<T>(vectorField, domain, 1);
    maskedComputeNorm(vectorField, *norm, mask, domain);
    return norm;
}


/* *************** Squared vector-norm of each cell in the field ******** */

template<typename T>
void computeNormSqr(MultiNTensorField2D<T>& vectorField,
                    MultiNTensorField2D<T>& normSqr, Box2D domain)
{
    applyProcessingFunctional (
            new ComputeNTensorNormSqrFunctional2D<T>, domain, normSqr, vectorField );
}

template<typename T>
void maskedComputeNormSqr(MultiNTensorField2D<T>& vectorField,
                          MultiNTensorField2D<T>& normSqr, MultiNTensorField2D<int>& mask, Box2D domain)
{
    applyProcessingFunctional (
            new MaskedComputeNTensorNormSqrFunctional2D<T>, domain, normSqr, vectorField, mask );
}

template<typename T>
MultiNTensorField2D<T>* computeNormSqr(MultiNTensorField2D<T>& vectorField, Box2D domain)
{
    MultiNTensorField2D<T>* normSqr
        = generateMultiNTensorField<T>(vectorField, domain, 1);
    computeNormSqr(vectorField, *normSqr, domain);
    return normSqr;
}

template<typename T>
MultiNTensorField2D<T>* maskedComputeNormSqr(MultiNTensorField2D<T>& vectorField, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<T>* normSqr
        = generateMultiNTensorField<T>(vectorField, domain, 1);
    maskedComputeNormSqr(vectorField, *normSqr, mask, domain);
    return normSqr;
}


/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template<typename T>
void computeSymmetricTensorNorm(MultiNTensorField2D<T>& tensorField,
                                MultiNTensorField2D<T>& norm,
                                Box2D domain)
{
    PLB_PRECONDITION( tensorField.getNdim()==3 );
    PLB_PRECONDITION( norm.getNdim()==1 );
    applyProcessingFunctional (
            new ComputeSymmetricNTensorNormFunctional2D<T>, domain, norm, tensorField );
}

template<typename T>
void maskedComputeSymmetricTensorNorm(MultiNTensorField2D<T>& tensorField,
                                      MultiNTensorField2D<T>& norm, MultiNTensorField2D<int>& mask,
                                      Box2D domain)
{
    PLB_PRECONDITION( tensorField.getNdim()==3 );
    PLB_PRECONDITION( norm.getNdim()==1 );
    applyProcessingFunctional (
            new MaskedComputeSymmetricNTensorNormFunctional2D<T>, domain, norm, tensorField, mask );
}

template<typename T>
MultiNTensorField2D<T>* computeSymmetricTensorNorm (
        MultiNTensorField2D<T>& tensorField, Box2D domain )
{
    MultiNTensorField2D<T>* norm
        = generateMultiNTensorField<T>(tensorField, domain, 1);
    computeSymmetricTensorNorm(tensorField, *norm, domain);
    return norm;
}

template<typename T>
MultiNTensorField2D<T>* maskedComputeSymmetricTensorNorm (
        MultiNTensorField2D<T>& tensorField, MultiNTensorField2D<int>& mask, Box2D domain )
{
    MultiNTensorField2D<T>* norm
        = generateMultiNTensorField<T>(tensorField, domain, 1);
    maskedComputeSymmetricTensorNorm(tensorField, *norm, mask, domain);
    return norm;
}


/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template<typename T>
void computeSymmetricTensorNormSqr(MultiNTensorField2D<T>& tensorField,
                                   MultiNTensorField2D<T>& normSqr, Box2D domain)
{
    PLB_PRECONDITION( tensorField.getNdim()==3 );
    PLB_PRECONDITION( normSqr.getNdim()==1 );
    applyProcessingFunctional (
            new ComputeSymmetricNTensorNormSqrFunctional2D<T>,
            domain, normSqr, tensorField );
}

template<typename T>
void maskedComputeSymmetricTensorNormSqr(MultiNTensorField2D<T>& tensorField,
                                         MultiNTensorField2D<T>& normSqr, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( tensorField.getNdim()==3 );
    PLB_PRECONDITION( normSqr.getNdim()==1 );
    applyProcessingFunctional (
            new MaskedComputeSymmetricNTensorNormSqrFunctional2D<T>,
            domain, normSqr, tensorField, mask );
}

template<typename T>
MultiNTensorField2D<T>* computeSymmetricTensorNormSqr (
        MultiNTensorField2D<T>& tensorField, Box2D domain )
{
    MultiNTensorField2D<T>* normSqr
        = generateMultiNTensorField<T>(tensorField, domain, 1);
    computeSymmetricTensorNormSqr(tensorField, *normSqr, domain);
    return normSqr;
}

template<typename T>
MultiNTensorField2D<T>* maskedComputeSymmetricTensorNormSqr (
        MultiNTensorField2D<T>& tensorField, MultiNTensorField2D<int>& mask, Box2D domain )
{
    MultiNTensorField2D<T>* normSqr
        = generateMultiNTensorField<T>(tensorField, domain, 1);
    maskedComputeSymmetricTensorNormSqr(tensorField, *normSqr, mask, domain);
    return normSqr;
}

/* *************** Trace of each symmetric tensor of a field ************ */

template<typename T>
void computeSymmetricTensorTrace(MultiNTensorField2D<T>& tensorField,
                                 MultiNTensorField2D<T>& trace,
                                 Box2D domain)
{
    PLB_PRECONDITION( tensorField.getNdim()==3 );
    PLB_PRECONDITION( trace.getNdim()==1 );
    applyProcessingFunctional (
            new ComputeSymmetricNTensorTraceFunctional2D<T>, domain, trace, tensorField );
}

template<typename T>
void maskedComputeSymmetricTensorTrace(MultiNTensorField2D<T>& tensorField,
                                       MultiNTensorField2D<T>& trace, MultiNTensorField2D<int>& mask,
                                       Box2D domain)
{
    PLB_PRECONDITION( tensorField.getNdim()==3 );
    PLB_PRECONDITION( trace.getNdim()==1 );
    applyProcessingFunctional (
            new MaskedComputeSymmetricNTensorTraceFunctional2D<T>, domain, trace, tensorField, mask );
}

template<typename T>
MultiNTensorField2D<T>* computeSymmetricTensorTrace(MultiNTensorField2D<T>& tensorField,
                                                    Box2D domain)
{
    MultiNTensorField2D<T>* trace
        = generateMultiNTensorField<T>(tensorField, domain, 1);
    computeSymmetricTensorTrace(tensorField, *trace, domain);
    return trace;
}

template<typename T>
MultiNTensorField2D<T>* maskedComputeSymmetricTensorTrace(MultiNTensorField2D<T>& tensorField, MultiNTensorField2D<int>& mask,
                                                          Box2D domain)
{
    MultiNTensorField2D<T>* trace
        = generateMultiNTensorField<T>(tensorField, domain, 1);
    maskedComputeSymmetricTensorTrace(tensorField, *trace, mask, domain);
    return trace;
}


/* *************** Vorticity from Velocity field *********************** */

template<typename T>
void computeVorticity( MultiNTensorField2D<T>& velocity,
                       MultiNTensorField2D<T>& vorticity, Box2D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxNTensorVorticityFunctional2D<T>,
            domain, vorticity, velocity, envelopeWidth );
}

template<typename T>
void maskedComputeVorticity( MultiNTensorField2D<T>& velocity,
                             MultiNTensorField2D<T>& vorticity, MultiNTensorField2D<int>& mask, Box2D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new MaskedBoxNTensorVorticityFunctional2D<T>,
            domain, vorticity, velocity, mask, envelopeWidth );
}

template<typename T>
MultiNTensorField2D<T>* computeVorticity(MultiNTensorField2D<T>& velocity, Box2D domain)
{
    MultiNTensorField2D<T>* vorticity
        = generateMultiNTensorField<T>(velocity, domain, 1);
    computeVorticity(velocity, *vorticity, domain);
    return vorticity;
}

template<typename T>
MultiNTensorField2D<T>* maskedComputeVorticity(MultiNTensorField2D<T>& velocity, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<T>* vorticity
        = generateMultiNTensorField<T>(velocity, domain, 1);
    maskedComputeVorticity(velocity, *vorticity, mask, domain);
    return vorticity;
}


/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkVorticity( MultiNTensorField2D<T>& velocity,
                           MultiNTensorField2D<T>& vorticity, Box2D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    applyProcessingFunctional (
            new BoxBulkNTensorVorticityFunctional2D<T>, domain, vorticity, velocity );
}

template<typename T>
void maskedComputeBulkVorticity( MultiNTensorField2D<T>& velocity,
                                 MultiNTensorField2D<T>& vorticity, MultiNTensorField2D<int>& mask, Box2D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    applyProcessingFunctional (
            new MaskedBoxBulkNTensorVorticityFunctional2D<T>, domain, vorticity, velocity, mask );
}

template<typename T>
MultiNTensorField2D<T>* computeBulkVorticity(MultiNTensorField2D<T>& velocity, Box2D domain)
{
    MultiNTensorField2D<T>* vorticity
        = generateMultiNTensorField<T>(velocity, domain, 1);
    computeBulkVorticity(velocity, *vorticity, domain);
    return vorticity;
}

template<typename T>
MultiNTensorField2D<T>* maskedComputeBulkVorticity(MultiNTensorField2D<T>& velocity, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<T>* vorticity
        = generateMultiNTensorField<T>(velocity, domain, 1);
    maskedComputeBulkVorticity(velocity, *vorticity, mask, domain);
    return vorticity;
}


/* *************** Strain Rate from Velocity field ********************* */

template<typename T>
void computeStrainRate( MultiNTensorField2D<T>& velocity,
                        MultiNTensorField2D<T>& S, Box2D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( S.getNdim()==3 );
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxNTensorStrainRateFunctional2D<T>, domain, velocity, S, envelopeWidth );
}

template<typename T>
void maskedComputeStrainRate( MultiNTensorField2D<T>& velocity,
                              MultiNTensorField2D<T>& S, MultiNTensorField2D<int>& mask, Box2D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( S.getNdim()==3 );
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new MaskedBoxNTensorStrainRateFunctional2D<T>, domain, velocity, S, mask, envelopeWidth );
}

template<typename T>
MultiNTensorField2D<T>* computeStrainRate(MultiNTensorField2D<T>& velocity, Box2D domain)
{
    MultiNTensorField2D<T>* S
        = generateMultiNTensorField<T>(velocity, domain, 3);
    computeStrainRate(velocity, *S, domain);
    return S;
}

template<typename T>
MultiNTensorField2D<T>* maskedComputeStrainRate(MultiNTensorField2D<T>& velocity, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<T>* S
        = generateMultiNTensorField<T>(velocity, domain, 3);
    maskedComputeStrainRate(velocity, *S, mask, domain);
    return S;
}


/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkStrainRate ( MultiNTensorField2D<T>& velocity,
                             MultiNTensorField2D<T>& S, Box2D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( S.getNdim()==3 );
    applyProcessingFunctional (
            new BoxBulkNTensorStrainRateFunctional2D<T>, domain, velocity, S );
}

template<typename T>
void maskedComputeBulkStrainRate ( MultiNTensorField2D<T>& velocity,
                                   MultiNTensorField2D<T>& S, MultiNTensorField2D<int>& mask, Box2D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( S.getNdim()==3 );
    applyProcessingFunctional (
            new MaskedBoxBulkNTensorStrainRateFunctional2D<T>, domain, velocity, S, mask );
}

template<typename T>
MultiNTensorField2D<T>* computeBulkStrainRate(MultiNTensorField2D<T>& velocity, Box2D domain)
{
    MultiNTensorField2D<T>* S
        = generateMultiNTensorField<T>(velocity, domain, 3);
    computeBulkStrainRate(velocity, *S, domain);
    return S;
}

template<typename T>
MultiNTensorField2D<T>* maskedComputeBulkStrainRate(MultiNTensorField2D<T>& velocity, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<T>* S
        = generateMultiNTensorField<T>(velocity, domain, 3);
    maskedComputeBulkStrainRate(velocity, *S, mask, domain);
    return S;
}


/* *************** MultiNTensorField - MultiNTensorField operations *************** */

template<typename T>
void add(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
         MultiNTensorField2D<T>& result, Box2D domain)
{
    std::vector<MultiNTensorField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_plus_B_NTensor2D<T>, domain, fields );
}

template<typename T>
void maskedAdd(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
         MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    std::vector<MultiNTensorField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Masked_A_plus_B_NTensor2D<T>, domain, fields, mask );
}

template<typename T>
MultiNTensorField2D<T>* add(MultiNTensorField2D<T>& A,
                            MultiNTensorField2D<T>& B, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<T>* result
        = generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    add(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<T>* maskedAdd(MultiNTensorField2D<T>& A,
                                  MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<T>* result
        = generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    maskedAdd(A, B, *result, mask, domain);
    return result;
}


template<typename T>
void add( MultiNTensorField2D<T>& A, T* alpha, int size,
          MultiNTensorField2D<T>& result, Box2D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_plus_alpha_NTensor2D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedAdd( MultiNTensorField2D<T>& A, T* alpha, int size,
                MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_plus_alpha_NTensor2D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField2D<T>* add( MultiNTensorField2D<T>& A,
                             T* alpha, int size, Box2D domain )
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    add(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<T>* maskedAdd( MultiNTensorField2D<T>& A,
                                   T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain )
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedAdd(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void subtract(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
              MultiNTensorField2D<T>& result, Box2D domain)
{
    std::vector<MultiNTensorField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_minus_B_NTensor2D<T>, domain, fields );
}

template<typename T>
void maskedSubtract(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                    MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    std::vector<MultiNTensorField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Masked_A_minus_B_NTensor2D<T>, domain, fields, mask);
}

template<typename T>
MultiNTensorField2D<T>* subtract( MultiNTensorField2D<T>& A,
                                  MultiNTensorField2D<T>& B, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<T>* result =
        generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    subtract(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<T>* maskedSubtract( MultiNTensorField2D<T>& A,
                                        MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<T>* result =
        generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    maskedSubtract(A, B, *result, mask, domain);
    return result;
}


template<typename T>
void subtract(MultiNTensorField2D<T>& A, T* alpha, int size,
              MultiNTensorField2D<T>& result, Box2D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_minus_alpha_NTensor2D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedSubtract(MultiNTensorField2D<T>& A, T* alpha, int size,
                    MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_minus_alpha_NTensor2D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField2D<T>* subtract( MultiNTensorField2D<T>& A,
                                  T* alpha, int size, Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    subtract(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<T>* maskedSubtract( MultiNTensorField2D<T>& A,
                                        T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedSubtract(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void subtract(T* alpha, int size, MultiNTensorField2D<T>& A,
              MultiNTensorField2D<T>& result, Box2D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Alpha_minus_A_NTensor2D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedSubtract(T* alpha, int size, MultiNTensorField2D<T>& A,
                    MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_Alpha_minus_A_NTensor2D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField2D<T>* subtract(T* alpha, int size,
                                 MultiNTensorField2D<T>& A, Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    subtract(alpha, size, A, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<T>* maskedSubtract(T* alpha, int size,
                                       MultiNTensorField2D<T>& A, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedSubtract(alpha, size, A, *result, mask, domain);
    return result;
}

template<typename T>
void multiply(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
              MultiNTensorField2D<T>& result, Box2D domain)
{
    std::vector<MultiNTensorField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_times_B_NTensor2D<T>, domain, fields );
}

template<typename T>
void maskedMultiply(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                    MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    std::vector<MultiNTensorField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Masked_A_times_B_NTensor2D<T>, domain, fields, mask );
}

template<typename T>
MultiNTensorField2D<T>* multiply(MultiNTensorField2D<T>& A,
                                 MultiNTensorField2D<T>& B, Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    multiply(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<T>* maskedMultiply(MultiNTensorField2D<T>& A,
                                       MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedMultiply(A, B, *result, mask, domain);
    return result;
}


template<typename T>
void multiply(MultiNTensorField2D<T>& A, T* alpha, int size,
              MultiNTensorField2D<T>& result, Box2D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_times_alpha_NTensor2D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedMultiply(MultiNTensorField2D<T>& A, T* alpha, int size,
                    MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_times_alpha_NTensor2D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField2D<T>* multiply(MultiNTensorField2D<T>& A,
                                 T* alpha, int size, Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    multiply(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<T>* maskedMultiply(MultiNTensorField2D<T>& A,
                                       T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedMultiply(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void divide(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
            MultiNTensorField2D<T>& result, Box2D domain)
{
    std::vector<MultiNTensorField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_dividedBy_B_NTensor2D<T>, domain, fields );
}

template<typename T>
void maskedDivide(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                  MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    std::vector<MultiNTensorField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Masked_A_dividedBy_B_NTensor2D<T>, domain, fields, mask );
}

template<typename T>
MultiNTensorField2D<T>* divide(MultiNTensorField2D<T>& A,
                               MultiNTensorField2D<T>& B, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<T>* result =
        generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    divide(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<T>* maskedDivide(MultiNTensorField2D<T>& A,
                                     MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<T>* result =
        generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    maskedDivide(A, B, *result, mask, domain);
    return result;
}


template<typename T>
void divide(MultiNTensorField2D<T>& A, T* alpha, int size,
            MultiNTensorField2D<T>& result, Box2D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_dividedBy_alpha_NTensor2D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedDivide(MultiNTensorField2D<T>& A, T* alpha, int size,
                  MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_dividedBy_alpha_NTensor2D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField2D<T>* divide(MultiNTensorField2D<T>& A,
                               T* alpha, int size, Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    divide(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<T>* maskedDivide(MultiNTensorField2D<T>& A,
                                     T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedDivide(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void divide(T* alpha, int size, MultiNTensorField2D<T>& A,
            MultiNTensorField2D<T>& result, Box2D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Alpha_dividedBy_A_NTensor2D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedDivide(T* alpha, int size, MultiNTensorField2D<T>& A,
                  MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_Alpha_dividedBy_A_NTensor2D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField2D<T>* divide(T* alpha, int size, MultiNTensorField2D<T>& A,
                               Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    divide(alpha, size, A, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<T>* maskedDivide(T* alpha, int size, MultiNTensorField2D<T>& A, MultiNTensorField2D<int>& mask,
                                     Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedDivide(alpha, size, A, *result, mask, domain);
    return result;
}


template<typename T>
void toThePower(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                MultiNTensorField2D<T>& result, Box2D domain)
{
    std::vector<MultiNTensorField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_toThePower_B_NTensor2D<T>, domain, fields );
}

template<typename T>
void maskedToThePower(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                      MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    std::vector<MultiNTensorField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Masked_A_toThePower_B_NTensor2D<T>, domain, fields, mask );
}

template<typename T>
MultiNTensorField2D<T>* toThePower(MultiNTensorField2D<T>& A,
                                   MultiNTensorField2D<T>& B, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<T>* result =
        generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    toThePower(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<T>* maskedToThePower(MultiNTensorField2D<T>& A,
                                         MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<T>* result =
        generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    maskedToThePower(A, B, *result, mask, domain);
    return result;
}


template<typename T>
void toThePower(MultiNTensorField2D<T>& A, T* alpha, int size,
                MultiNTensorField2D<T>& result, Box2D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_toThePower_alpha_NTensor2D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedToThePower(MultiNTensorField2D<T>& A, T* alpha, int size,
                      MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_toThePower_alpha_NTensor2D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField2D<T>* toThePower(MultiNTensorField2D<T>& A,
                                   T* alpha, int size, Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    toThePower(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<T>* maskedToThePower(MultiNTensorField2D<T>& A,
                                         T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedToThePower(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void toThePower(T* alpha, int size, MultiNTensorField2D<T>& A,
                MultiNTensorField2D<T>& result, Box2D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Alpha_toThePower_A_NTensor2D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedToThePower(T* alpha, int size, MultiNTensorField2D<T>& A,
                      MultiNTensorField2D<T>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_Alpha_toThePower_A_NTensor2D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField2D<T>* toThePower(T* alpha, int size, MultiNTensorField2D<T>& A,
                                   Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    toThePower(alpha, size, A, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<T>* maskedToThePower(T* alpha, int size, MultiNTensorField2D<T>& A, MultiNTensorField2D<int>& mask,
                                         Box2D domain)
{
    MultiNTensorField2D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedToThePower(alpha, size, A, *result, mask, domain);
    return result;
}


template<typename T>
void equals(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
            MultiNTensorField2D<int>& result, Box2D domain)
{
    std::vector<MultiBlock2D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_equals_B_NTensor2D<T>, domain, fields );
}

template<typename T>
void maskedEquals(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                  MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    std::vector<MultiBlock2D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_A_equals_B_NTensor2D<T>, domain, fields );
}

template<typename T>
MultiNTensorField2D<int>* equals(MultiNTensorField2D<T>& A,
                                 MultiNTensorField2D<T>& B, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    equals(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<int>* maskedEquals(MultiNTensorField2D<T>& A,
                                       MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    maskedEquals(A, B, *result, mask, domain);
    return result;
}

template<typename T>
void equals( MultiNTensorField2D<T>& A, T* alpha, int size,
             MultiNTensorField2D<int>& result, Box2D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_equals_alpha_NTensor2D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedEquals( MultiNTensorField2D<T>& A, T* alpha, int size,
                   MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_equals_alpha_NTensor2D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField2D<int>* equals(MultiNTensorField2D<T>& A,
                                 T* alpha, int size, Box2D domain)
{
    MultiNTensorField2D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    equals(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<int>* maskedEquals(MultiNTensorField2D<T>& A,
                                       T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    maskedEquals(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void lessThan(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
              MultiNTensorField2D<int>& result, Box2D domain)
{
    std::vector<MultiBlock2D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_lessThan_B_NTensor2D<T>, domain, fields );
}

template<typename T>
void maskedLessThan(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                    MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    std::vector<MultiBlock2D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_A_lessThan_B_NTensor2D<T>, domain, fields );
}

template<typename T>
MultiNTensorField2D<int>* lessThan(MultiNTensorField2D<T>& A,
                                   MultiNTensorField2D<T>& B, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    lessThan(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<int>* maskedLessThan(MultiNTensorField2D<T>& A,
                                         MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    maskedLessThan(A, B, *result, mask, domain);
    return result;
}

template<typename T>
void lessThan( MultiNTensorField2D<T>& A, T* alpha, int size,
               MultiNTensorField2D<int>& result, Box2D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_lessThan_alpha_NTensor2D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedLessThan( MultiNTensorField2D<T>& A, T* alpha, int size,
                     MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_lessThan_alpha_NTensor2D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField2D<int>* lessThan(MultiNTensorField2D<T>& A,
                                   T* alpha, int size, Box2D domain)
{
    MultiNTensorField2D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    lessThan(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<int>* maskedLessThan(MultiNTensorField2D<T>& A,
                                         T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    maskedLessThan(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void lessEqual(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
               MultiNTensorField2D<int>& result, Box2D domain)
{
    std::vector<MultiBlock2D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_lessEqual_B_NTensor2D<T>, domain, fields );
}

template<typename T>
void maskedLessEqual(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                     MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    std::vector<MultiBlock2D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_A_lessEqual_B_NTensor2D<T>, domain, fields );
}

template<typename T>
MultiNTensorField2D<int>* lessEqual(MultiNTensorField2D<T>& A,
                                    MultiNTensorField2D<T>& B, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    lessEqual(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<int>* maskedLessEqual(MultiNTensorField2D<T>& A,
                                          MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    maskedLessEqual(A, B, *result, mask, domain);
    return result;
}

template<typename T>
void lessEqual( MultiNTensorField2D<T>& A, T* alpha, int size,
                MultiNTensorField2D<int>& result, Box2D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_lessEqual_alpha_NTensor2D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedLessEqual( MultiNTensorField2D<T>& A, T* alpha, int size,
                      MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_lessEqual_alpha_NTensor2D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField2D<int>* lessEqual(MultiNTensorField2D<T>& A,
                                    T* alpha, int size, Box2D domain)
{
    MultiNTensorField2D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    lessEqual(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<int>* maskedLessEqual(MultiNTensorField2D<T>& A,
                                          T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    maskedLessEqual(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void greaterThan(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                 MultiNTensorField2D<int>& result, Box2D domain)
{
    std::vector<MultiBlock2D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_greaterThan_B_NTensor2D<T>, domain, fields );
}

template<typename T>
void maskedGreaterThan(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                       MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    std::vector<MultiBlock2D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_A_greaterThan_B_NTensor2D<T>, domain, fields );
}

template<typename T>
MultiNTensorField2D<int>* greaterThan(MultiNTensorField2D<T>& A,
                                      MultiNTensorField2D<T>& B, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    greaterThan(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<int>* maskedGreaterThan(MultiNTensorField2D<T>& A,
                                            MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    maskedGreaterThan(A, B, *result, mask, domain);
    return result;
}

template<typename T>
void greaterThan( MultiNTensorField2D<T>& A, T* alpha, int size,
                  MultiNTensorField2D<int>& result, Box2D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_greaterThan_alpha_NTensor2D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedGreaterThan( MultiNTensorField2D<T>& A, T* alpha, int size,
                        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_greaterThan_alpha_NTensor2D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField2D<int>* greaterThan(MultiNTensorField2D<T>& A,
                                      T* alpha, int size, Box2D domain)
{
    MultiNTensorField2D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    greaterThan(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<int>* maskedGreaterThan(MultiNTensorField2D<T>& A,
                                      T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    maskedGreaterThan(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void greaterEqual(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                  MultiNTensorField2D<int>& result, Box2D domain)
{
    std::vector<MultiBlock2D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_greaterEqual_B_NTensor2D<T>, domain, fields );
}

template<typename T>
void maskedGreaterEqual(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                        MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    std::vector<MultiBlock2D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_A_greaterEqual_B_NTensor2D<T>, domain, fields );
}

template<typename T>
MultiNTensorField2D<int>* greaterEqual(MultiNTensorField2D<T>& A,
                                       MultiNTensorField2D<T>& B, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    greaterEqual(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<int>* maskedGreaterEqual(MultiNTensorField2D<T>& A,
                                             MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    maskedGreaterEqual(A, B, *result, mask, domain);
    return result;
}

template<typename T>
void greaterEqual( MultiNTensorField2D<T>& A, T* alpha, int size,
                   MultiNTensorField2D<int>& result, Box2D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_greaterEqual_alpha_NTensor2D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedGreaterEqual( MultiNTensorField2D<T>& A, T* alpha, int size,
                         MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_greaterEqual_alpha_NTensor2D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField2D<int>* greaterEqual(MultiNTensorField2D<T>& A,
                                       T* alpha, int size, Box2D domain)
{
    MultiNTensorField2D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    greaterEqual(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<int>* maskedGreaterEqual(MultiNTensorField2D<T>& A,
                                             T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    maskedGreaterEqual(A, alpha, size, *result, mask, domain);
    return result;
}

template<typename T>
void logicalAnd(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                MultiNTensorField2D<int>& result, Box2D domain)
{
    std::vector<MultiBlock2D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_and_B_NTensor2D<T>, domain, fields );
}

template<typename T>
void maskedLogicalAnd(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                      MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    std::vector<MultiBlock2D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_A_and_B_NTensor2D<T>, domain, fields );
}

template<typename T>
MultiNTensorField2D<int>* logicalAnd(MultiNTensorField2D<T>& A,
                                     MultiNTensorField2D<T>& B, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    logicalAnd(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<int>* maskedLogicalAnd(MultiNTensorField2D<T>& A,
                                           MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    maskedLogicalAnd(A, B, *result, mask, domain);
    return result;
}


template<typename T>
void logicalOr(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
               MultiNTensorField2D<int>& result, Box2D domain)
{
    std::vector<MultiBlock2D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_or_B_NTensor2D<T>, domain, fields );
}

template<typename T>
void maskedLogicalOr(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B,
                     MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain)
{
    std::vector<MultiBlock2D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_A_or_B_NTensor2D<T>, domain, fields );
}

template<typename T>
MultiNTensorField2D<int>* logicalOr(MultiNTensorField2D<T>& A,
                                    MultiNTensorField2D<T>& B, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    logicalOr(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<int>* maskedLogicalOr(MultiNTensorField2D<T>& A,
                                          MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField2D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    maskedLogicalOr(A, B, *result, mask, domain);
    return result;
}

template<typename T>
void negate( MultiNTensorField2D<T>& A,
             MultiNTensorField2D<int>& result, Box2D domain )
{
    applyProcessingFunctional (
            new Not_A_NTensor2D<T>(), domain, A, result );
}

template<typename T>
void maskedNegate( MultiNTensorField2D<T>& A,
                   MultiNTensorField2D<int>& result, MultiNTensorField2D<int>& mask, Box2D domain )
{
    applyProcessingFunctional (
            new Masked_Not_A_NTensor2D<T>(), domain, A, result, mask );
}

template<typename T>
MultiNTensorField2D<int>* negate(MultiNTensorField2D<T>& A, Box2D domain)
{
    MultiNTensorField2D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    negate(A, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField2D<int>* maskedNegate(MultiNTensorField2D<T>& A, MultiNTensorField2D<int>& mask, Box2D domain)
{
    MultiNTensorField2D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    maskedNegate(A, *result, mask, domain);
    return result;
}


/* *************** NTensorField - NTensorField inplace operations *************** */

template<typename T>
void assignInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, Box2D domain) {
    applyProcessingFunctional (
            new AssignNTensorFunctional2D<T>, domain, A, B );
}

template<typename T>
void maskedAssignInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain) {
    applyProcessingFunctional (
            new MaskedAssignNTensorFunctional2D<T>, domain, A, B, mask );
}


template<typename T>
void assignInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, Box2D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new IniConstNTensorFunctional2D<T>(alphaVect), domain, A );
}

template<typename T>
void maskedAssignInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new MaskedIniConstNTensorFunctional2D<T>(alphaVect), domain, A, mask );
}


template<typename T>
void addInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, Box2D domain) {
    applyProcessingFunctional (
            new A_plus_B_inplace_NTensor2D<T>, domain, A, B );
}

template<typename T>
void maskedAddInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain) {
    applyProcessingFunctional (
            new Masked_A_plus_B_inplace_NTensor2D<T>, domain, A, B, mask );
}


template<typename T>
void addInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, Box2D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_plus_alpha_inplace_NTensor2D<T>(alphaVect), domain, A );
}

template<typename T>
void maskedAddInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_plus_alpha_inplace_NTensor2D<T>(alphaVect), domain, A, mask );
}


template<typename T>
void subtractInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, Box2D domain) {
    applyProcessingFunctional (
            new A_minus_B_inplace_NTensor2D<T>, domain, A, B );
}

template<typename T>
void maskedSubtractInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain) {
    applyProcessingFunctional (
            new Masked_A_minus_B_inplace_NTensor2D<T>, domain, A, B, mask );
}


template<typename T>
void subtractInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, Box2D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_minus_alpha_inplace_NTensor2D<T>(alphaVect), domain, A );
}

template<typename T>
void maskedSubtractInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_minus_alpha_inplace_NTensor2D<T>(alphaVect), domain, A, mask );
}


template<typename T>
void multiplyInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, Box2D domain) {
    applyProcessingFunctional (
            new A_times_B_inplace_NTensor2D<T>, domain, A, B );
}

template<typename T>
void maskedMultiplyInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain) {
    applyProcessingFunctional (
            new Masked_A_times_B_inplace_NTensor2D<T>, domain, A, B, mask );
}


template<typename T>
void multiplyInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, Box2D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_times_alpha_inplace_NTensor2D<T>(alphaVect), domain, A );
}

template<typename T>
void maskedMultiplyInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_times_alpha_inplace_NTensor2D<T>(alphaVect), domain, A, mask );
}


template<typename T>
void divideInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, Box2D domain) {
    applyProcessingFunctional (
            new A_dividedBy_B_inplace_NTensor2D<T>, domain, A, B );
}

template<typename T>
void maskedDivideInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain) {
    applyProcessingFunctional (
            new Masked_A_dividedBy_B_inplace_NTensor2D<T>, domain, A, B, mask );
}


template<typename T>
void divideInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, Box2D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_dividedBy_alpha_inplace_NTensor2D<T>(alphaVect), domain, A );
}

template<typename T>
void maskedDivideInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_dividedBy_alpha_inplace_NTensor2D<T>(alphaVect), domain, A, mask );
}


template<typename T>
void toThePowerInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, Box2D domain) {
    applyProcessingFunctional (
            new A_toThePower_B_inplace_NTensor2D<T>, domain, A, B );
}

template<typename T>
void maskedToThePowerInPlace(MultiNTensorField2D<T>& A, MultiNTensorField2D<T>& B, MultiNTensorField2D<int>& mask, Box2D domain) {
    applyProcessingFunctional (
            new Masked_A_toThePower_B_inplace_NTensor2D<T>, domain, A, B, mask );
}


template<typename T>
void toThePowerInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, Box2D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_toThePower_alpha_inplace_NTensor2D<T>(alphaVect), domain, A );
}

template<typename T>
void maskedToThePowerInPlace(MultiNTensorField2D<T>& A, T* alpha, int size, MultiNTensorField2D<int>& mask, Box2D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_toThePower_alpha_inplace_NTensor2D<T>(alphaVect), domain, A, mask );
}

}  // namespace plb

#endif  // SWIG_DATA_ANALYSIS_WRAPPER_2D_HH
