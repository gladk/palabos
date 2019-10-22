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

#ifndef SWIG_DATA_ANALYSIS_WRAPPER_3D_HH
#define SWIG_DATA_ANALYSIS_WRAPPER_3D_HH

#include "plbWrapper/block/dataAnalysisWrapper3D.h"
#include "plbWrapper/block/dataAnalysisFunctional3D.h"
#include "plbWrapper/block/dataInitializerFunctional3D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "multiBlock/multiBlockGenerator3D.h"

namespace plb {


/* *************** Reductive functions ******************************* */

template<typename T>
void computeAverage(MultiNTensorField3D<T>& vectorField, Box3D domain, T* result, int size)
{
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    BoxNTensorSumFunctional3D<T> functional(vectorField.getNdim());
    applyProcessingFunctional(functional, domain, vectorField);
    std::vector<T> const& sumVector = functional.getSumVector();
    for (plint iDim=0; iDim<vectorField.getNdim(); ++iDim) {
        result[iDim] = sumVector[iDim] / (T) domain.nCells();
    }
}

template<typename T>
void maskedComputeAverage(MultiNTensorField3D<T>& vectorField, MultiNTensorField3D<int>& mask, Box3D domain, T* result, int size)
{
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    MaskedBoxNTensorSumFunctional3D<T> functional(vectorField.getNdim());
    applyProcessingFunctional(functional, domain, vectorField, mask);
    std::vector<T> const& sumVector = functional.getSumVector();
    for (plint iDim=0; iDim<vectorField.getNdim(); ++iDim) {
        result[iDim] = sumVector[iDim] / (T) domain.nCells();
    }
}

template<typename T>
void computeMin(MultiNTensorField3D<T>& vectorField, Box3D domain, T* result, int size)
{
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    BoxNTensorMinFunctional3D<T> functional(vectorField.getNdim());
    applyProcessingFunctional(functional, domain, vectorField);
    std::vector<T> const& minVector = functional.getMinVector();
    for (plint iDim=0; iDim<vectorField.getNdim(); ++iDim) {
        result[iDim] = minVector[iDim];
    }
}

template<typename T>
void maskedComputeMin(MultiNTensorField3D<T>& vectorField, MultiNTensorField3D<int>& mask, Box3D domain, T* result, int size)
{
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    MaskedBoxNTensorMinFunctional3D<T> functional(vectorField.getNdim());
    applyProcessingFunctional(functional, domain, vectorField, mask);
    std::vector<T> const& minVector = functional.getMinVector();
    for (plint iDim=0; iDim<vectorField.getNdim(); ++iDim) {
        result[iDim] = minVector[iDim];
    }
}


template<typename T>
void computeMax(MultiNTensorField3D<T>& vectorField, Box3D domain, T* result, int size)
{
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    BoxNTensorMaxFunctional3D<T> functional(vectorField.getNdim());
    applyProcessingFunctional(functional, domain, vectorField);
    std::vector<T> const& maxVector = functional.getMaxVector();
    for (plint iDim=0; iDim<vectorField.getNdim(); ++iDim) {
        result[iDim] = maxVector[iDim];
    }
}

template<typename T>
void maskedComputeMax(MultiNTensorField3D<T>& vectorField, MultiNTensorField3D<int>& mask, Box3D domain, T* result, int size)
{
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    MaskedBoxNTensorMaxFunctional3D<T> functional(vectorField.getNdim());
    applyProcessingFunctional(functional, domain, vectorField, mask);
    std::vector<T> const& maxVector = functional.getMaxVector();
    for (plint iDim=0; iDim<vectorField.getNdim(); ++iDim) {
        result[iDim] = maxVector[iDim];
    }
}

template<typename T>
void computeBoundedAverage(MultiNTensorField3D<T>& vectorField, Box3D domain, T* result, int size) {
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    BoundedBoxNTensorSumFunctional3D<T> functional(vectorField.getNdim());
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
        MultiNTensorField3D<T>& vectorField, MultiNTensorField3D<int>& mask,
        Box3D domain, T* result, int size)
{
    PLB_PRECONDITION( size == (int)vectorField.getNdim() );
    BoundedMaskedBoxNTensorSumFunctional3D<T> functional(vectorField.getNdim());
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
void maskedCopy( MultiNTensorField3D<T1>& field,
                 MultiNTensorField3D<T2>& convertedField, MultiNTensorField3D<int>& mask,
                 Box3D domain)
{
    applyProcessingFunctional (
            new MaskedCopyConvertNTensorFunctional3D<T1,T2>, domain, field, convertedField, mask );
}


template<typename T1, typename T2>
MultiNTensorField3D<T2>* maskedCopyConvert( MultiNTensorField3D<T1>& field, MultiNTensorField3D<int>& mask,
                                            Box3D domain)
{
    MultiNTensorField3D<T2>* convertedField
        = generateMultiNTensorField<T2>(field, domain, field.getNdim());
    maskedCopy(field, *convertedField, mask, domain);
    return convertedField;
}

/* *************** Component (scalar-field) out of a tensor-field ****** */

template<typename T>
void extractComponent(MultiNTensorField3D<T>& vectorField, MultiNTensorField3D<T>& component,
                      Box3D domain, int iComponent)
{
    applyProcessingFunctional (
            new ExtractNTensorComponentFunctional3D<T>(iComponent),
            domain, component, vectorField );
}

template<typename T>
void maskedExtractComponent(MultiNTensorField3D<T>& vectorField, MultiNTensorField3D<T>& component, MultiNTensorField3D<int>& mask,
                            Box3D domain, int iComponent)
{
    applyProcessingFunctional (
            new MaskedExtractNTensorComponentFunctional3D<T>(iComponent),
            domain, component, vectorField, mask );
}

template<typename T>
MultiNTensorField3D<T>* extractComponent (
              MultiNTensorField3D<T>& vectorField,
              Box3D domain, int iComponent )
{
    MultiNTensorField3D<T>* component
        = generateMultiNTensorField<T>(vectorField, domain, 1);
    extractComponent(vectorField, *component, domain, iComponent);
    return component;
}

template<typename T>
MultiNTensorField3D<T>* maskedExtractComponent (
        MultiNTensorField3D<T>& vectorField, MultiNTensorField3D<int>& mask,
        Box3D domain, int iComponent )
{
    MultiNTensorField3D<T>* component
        = generateMultiNTensorField<T>(vectorField, domain, 1);
    maskedExtractComponent(vectorField, *component, mask, domain, iComponent);
    return component;
}

/* *************** Vector-norm of each cell in the field *************** */

template<typename T>
void computeNorm(MultiNTensorField3D<T>& vectorField,
                 MultiNTensorField3D<T>& norm, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeNTensorNormFunctional3D<T>, domain, norm, vectorField );
}

template<typename T>
void maskedComputeNorm(MultiNTensorField3D<T>& vectorField,
                 MultiNTensorField3D<T>& norm, MultiNTensorField3D<int>& mask, Box3D domain)
{
    applyProcessingFunctional (
            new MaskedComputeNTensorNormFunctional3D<T>, domain, norm, vectorField, mask );
}

template<typename T>
MultiNTensorField3D<T>* computeNorm(MultiNTensorField3D<T>& vectorField, Box3D domain)
{
    MultiNTensorField3D<T>* norm
        = generateMultiNTensorField<T>(vectorField, domain, 1);
    computeNorm(vectorField, *norm, domain);
    return norm;
}

template<typename T>
MultiNTensorField3D<T>* maskedComputeNorm(MultiNTensorField3D<T>& vectorField, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<T>* norm
        = generateMultiNTensorField<T>(vectorField, domain, 1);
    maskedComputeNorm(vectorField, *norm, mask, domain);
    return norm;
}


/* *************** Squared vector-norm of each cell in the field ******** */

template<typename T>
void computeNormSqr(MultiNTensorField3D<T>& vectorField,
                    MultiNTensorField3D<T>& normSqr, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeNTensorNormSqrFunctional3D<T>, domain, normSqr, vectorField );
}

template<typename T>
void maskedComputeNormSqr(MultiNTensorField3D<T>& vectorField,
                          MultiNTensorField3D<T>& normSqr, MultiNTensorField3D<int>& mask, Box3D domain)
{
    applyProcessingFunctional (
            new MaskedComputeNTensorNormSqrFunctional3D<T>, domain, normSqr, vectorField, mask );
}

template<typename T>
MultiNTensorField3D<T>* computeNormSqr(MultiNTensorField3D<T>& vectorField, Box3D domain)
{
    MultiNTensorField3D<T>* normSqr
        = generateMultiNTensorField<T>(vectorField, domain, 1);
    computeNormSqr(vectorField, *normSqr, domain);
    return normSqr;
}

template<typename T>
MultiNTensorField3D<T>* maskedComputeNormSqr(MultiNTensorField3D<T>& vectorField, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<T>* normSqr
        = generateMultiNTensorField<T>(vectorField, domain, 1);
    maskedComputeNormSqr(vectorField, *normSqr, mask, domain);
    return normSqr;
}


/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template<typename T>
void computeSymmetricTensorNorm(MultiNTensorField3D<T>& tensorField,
                                MultiNTensorField3D<T>& norm,
                                Box3D domain)
{
    PLB_PRECONDITION( tensorField.getNdim()==6 );
    PLB_PRECONDITION( norm.getNdim()==1 );
    applyProcessingFunctional (
            new ComputeSymmetricNTensorNormFunctional3D<T>, domain, norm, tensorField );
}

template<typename T>
void maskedComputeSymmetricTensorNorm(MultiNTensorField3D<T>& tensorField,
                                      MultiNTensorField3D<T>& norm, MultiNTensorField3D<int>& mask,
                                      Box3D domain)
{
    PLB_PRECONDITION( tensorField.getNdim()==6 );
    PLB_PRECONDITION( norm.getNdim()==1 );
    applyProcessingFunctional (
            new MaskedComputeSymmetricNTensorNormFunctional3D<T>, domain, norm, tensorField, mask );
}

template<typename T>
MultiNTensorField3D<T>* computeSymmetricTensorNorm (
        MultiNTensorField3D<T>& tensorField, Box3D domain )
{
    MultiNTensorField3D<T>* norm
        = generateMultiNTensorField<T>(tensorField, domain, 1);
    computeSymmetricTensorNorm(tensorField, *norm, domain);
    return norm;
}

template<typename T>
MultiNTensorField3D<T>* maskedComputeSymmetricTensorNorm (
        MultiNTensorField3D<T>& tensorField, MultiNTensorField3D<int>& mask, Box3D domain )
{
    MultiNTensorField3D<T>* norm
        = generateMultiNTensorField<T>(tensorField, domain, 1);
    maskedComputeSymmetricTensorNorm(tensorField, *norm, mask, domain);
    return norm;
}


/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template<typename T>
void computeSymmetricTensorNormSqr(MultiNTensorField3D<T>& tensorField,
                                   MultiNTensorField3D<T>& normSqr, Box3D domain)
{
    PLB_PRECONDITION( tensorField.getNdim()==6 );
    PLB_PRECONDITION( normSqr.getNdim()==1 );
    applyProcessingFunctional (
            new ComputeSymmetricNTensorNormSqrFunctional3D<T>,
            domain, normSqr, tensorField );
}

template<typename T>
void maskedComputeSymmetricTensorNormSqr(MultiNTensorField3D<T>& tensorField,
                                         MultiNTensorField3D<T>& normSqr, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( tensorField.getNdim()==6 );
    PLB_PRECONDITION( normSqr.getNdim()==1 );
    applyProcessingFunctional (
            new MaskedComputeSymmetricNTensorNormSqrFunctional3D<T>,
            domain, normSqr, tensorField, mask );
}

template<typename T>
MultiNTensorField3D<T>* computeSymmetricTensorNormSqr (
        MultiNTensorField3D<T>& tensorField, Box3D domain )
{
    MultiNTensorField3D<T>* normSqr
        = generateMultiNTensorField<T>(tensorField, domain, 1);
    computeSymmetricTensorNormSqr(tensorField, *normSqr, domain);
    return normSqr;
}

template<typename T>
MultiNTensorField3D<T>* maskedComputeSymmetricTensorNormSqr (
        MultiNTensorField3D<T>& tensorField, MultiNTensorField3D<int>& mask, Box3D domain )
{
    MultiNTensorField3D<T>* normSqr
        = generateMultiNTensorField<T>(tensorField, domain, 1);
    maskedComputeSymmetricTensorNormSqr(tensorField, *normSqr, mask, domain);
    return normSqr;
}

/* *************** Trace of each symmetric tensor of a field ************ */

template<typename T>
void computeSymmetricTensorTrace(MultiNTensorField3D<T>& tensorField,
                                 MultiNTensorField3D<T>& trace,
                                 Box3D domain)
{
    PLB_PRECONDITION( tensorField.getNdim()==6 );
    PLB_PRECONDITION( trace.getNdim()==1 );
    applyProcessingFunctional (
            new ComputeSymmetricNTensorTraceFunctional3D<T>, domain, trace, tensorField );
}

template<typename T>
void maskedComputeSymmetricTensorTrace(MultiNTensorField3D<T>& tensorField,
                                       MultiNTensorField3D<T>& trace, MultiNTensorField3D<int>& mask,
                                       Box3D domain)
{
    PLB_PRECONDITION( tensorField.getNdim()==6 );
    PLB_PRECONDITION( trace.getNdim()==1 );
    applyProcessingFunctional (
            new MaskedComputeSymmetricNTensorTraceFunctional3D<T>, domain, trace, tensorField, mask );
}

template<typename T>
MultiNTensorField3D<T>* computeSymmetricTensorTrace(MultiNTensorField3D<T>& tensorField,
                                                    Box3D domain)
{
    MultiNTensorField3D<T>* trace
        = generateMultiNTensorField<T>(tensorField, domain, 1);
    computeSymmetricTensorTrace(tensorField, *trace, domain);
    return trace;
}

template<typename T>
MultiNTensorField3D<T>* maskedComputeSymmetricTensorTrace(MultiNTensorField3D<T>& tensorField, MultiNTensorField3D<int>& mask,
                                                          Box3D domain)
{
    MultiNTensorField3D<T>* trace
        = generateMultiNTensorField<T>(tensorField, domain, 1);
    maskedComputeSymmetricTensorTrace(tensorField, *trace, mask, domain);
    return trace;
}


/* *************** Vorticity from Velocity field *********************** */

template<typename T>
void computeVorticity( MultiNTensorField3D<T>& velocity,
                       MultiNTensorField3D<T>& vorticity, Box3D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxNTensorVorticityFunctional3D<T>,
            domain, vorticity, velocity, envelopeWidth );
}

template<typename T>
void maskedComputeVorticity( MultiNTensorField3D<T>& velocity,
                             MultiNTensorField3D<T>& vorticity, MultiNTensorField3D<int>& mask, Box3D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new MaskedBoxNTensorVorticityFunctional3D<T>,
            domain, vorticity, velocity, mask, envelopeWidth );
}

template<typename T>
MultiNTensorField3D<T>* computeVorticity(MultiNTensorField3D<T>& velocity, Box3D domain)
{
    MultiNTensorField3D<T>* vorticity
        = generateMultiNTensorField<T>(velocity, domain, velocity.getNdim());
    computeVorticity(velocity, *vorticity, domain);
    return vorticity;
}

template<typename T>
MultiNTensorField3D<T>* maskedComputeVorticity(MultiNTensorField3D<T>& velocity, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<T>* vorticity
        = generateMultiNTensorField<T>(velocity, domain, velocity.getNdim());
    maskedComputeVorticity(velocity, *vorticity, mask, domain);
    return vorticity;
}


/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkVorticity( MultiNTensorField3D<T>& velocity,
                           MultiNTensorField3D<T>& vorticity, Box3D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    applyProcessingFunctional (
            new BoxBulkNTensorVorticityFunctional3D<T>, domain, vorticity, velocity );
}

template<typename T>
void maskedComputeBulkVorticity( MultiNTensorField3D<T>& velocity,
                                 MultiNTensorField3D<T>& vorticity, MultiNTensorField3D<int>& mask, Box3D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    applyProcessingFunctional (
            new MaskedBoxBulkNTensorVorticityFunctional3D<T>, domain, vorticity, velocity, mask );
}

template<typename T>
MultiNTensorField3D<T>* computeBulkVorticity(MultiNTensorField3D<T>& velocity, Box3D domain)
{
    MultiNTensorField3D<T>* vorticity
        = generateMultiNTensorField<T>(velocity, domain, velocity.getNdim());
    computeBulkVorticity(velocity, *vorticity, domain);
    return vorticity;
}

template<typename T>
MultiNTensorField3D<T>* maskedComputeBulkVorticity(MultiNTensorField3D<T>& velocity, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<T>* vorticity
        = generateMultiNTensorField<T>(velocity, domain, velocity.getNdim());
    maskedComputeBulkVorticity(velocity, *vorticity, mask, domain);
    return vorticity;
}


/* *************** Strain Rate from Velocity field ********************* */

template<typename T>
void computeStrainRate( MultiNTensorField3D<T>& velocity,
                        MultiNTensorField3D<T>& S, Box3D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( S.getNdim()==6 );
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxNTensorStrainRateFunctional3D<T>, domain, velocity, S, envelopeWidth );
}

template<typename T>
void maskedComputeStrainRate( MultiNTensorField3D<T>& velocity,
                              MultiNTensorField3D<T>& S, MultiNTensorField3D<int>& mask, Box3D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( S.getNdim()==6 );
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new MaskedBoxNTensorStrainRateFunctional3D<T>, domain, velocity, S, mask, envelopeWidth );
}

template<typename T>
MultiNTensorField3D<T>* computeStrainRate(MultiNTensorField3D<T>& velocity, Box3D domain)
{
    MultiNTensorField3D<T>* S
        = generateMultiNTensorField<T>(velocity, domain, 6);
    computeStrainRate(velocity, *S, domain);
    return S;
}

template<typename T>
MultiNTensorField3D<T>* maskedComputeStrainRate(MultiNTensorField3D<T>& velocity, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<T>* S
        = generateMultiNTensorField<T>(velocity, domain, 6);
    maskedComputeStrainRate(velocity, *S, mask, domain);
    return S;
}


/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkStrainRate ( MultiNTensorField3D<T>& velocity,
                             MultiNTensorField3D<T>& S, Box3D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( S.getNdim()==6 );
    applyProcessingFunctional (
            new BoxBulkNTensorStrainRateFunctional3D<T>, domain, velocity, S );
}

template<typename T>
void maskedComputeBulkStrainRate ( MultiNTensorField3D<T>& velocity,
                                   MultiNTensorField3D<T>& S, MultiNTensorField3D<int>& mask, Box3D domain )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( S.getNdim()==6 );
    applyProcessingFunctional (
            new MaskedBoxBulkNTensorStrainRateFunctional3D<T>, domain, velocity, S, mask );
}

template<typename T>
MultiNTensorField3D<T>* computeBulkStrainRate(MultiNTensorField3D<T>& velocity, Box3D domain)
{
    MultiNTensorField3D<T>* S
        = generateMultiNTensorField<T>(velocity, domain, 6);
    computeBulkStrainRate(velocity, *S, domain);
    return S;
}

template<typename T>
MultiNTensorField3D<T>* maskedComputeBulkStrainRate(MultiNTensorField3D<T>& velocity, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<T>* S
        = generateMultiNTensorField<T>(velocity, domain, 6);
    maskedComputeBulkStrainRate(velocity, *S, mask, domain);
    return S;
}


/* *************** MultiNTensorField - MultiNTensorField operations *************** */

template<typename T>
void add(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
         MultiNTensorField3D<T>& result, Box3D domain)
{
    std::vector<MultiNTensorField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_plus_B_NTensor3D<T>, domain, fields );
}

template<typename T>
void maskedAdd(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
         MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    std::vector<MultiNTensorField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Masked_A_plus_B_NTensor3D<T>, domain, fields, mask );
}

template<typename T>
MultiNTensorField3D<T>* add(MultiNTensorField3D<T>& A,
                            MultiNTensorField3D<T>& B, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<T>* result
        = generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    add(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<T>* maskedAdd(MultiNTensorField3D<T>& A,
                                  MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<T>* result
        = generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    maskedAdd(A, B, *result, mask, domain);
    return result;
}


template<typename T>
void add( MultiNTensorField3D<T>& A, T* alpha, int size,
          MultiNTensorField3D<T>& result, Box3D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_plus_alpha_NTensor3D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedAdd( MultiNTensorField3D<T>& A, T* alpha, int size,
                MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_plus_alpha_NTensor3D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField3D<T>* add( MultiNTensorField3D<T>& A,
                             T* alpha, int size, Box3D domain )
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    add(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<T>* maskedAdd( MultiNTensorField3D<T>& A,
                                   T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain )
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedAdd(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void subtract(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
              MultiNTensorField3D<T>& result, Box3D domain)
{
    std::vector<MultiNTensorField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_minus_B_NTensor3D<T>, domain, fields );
}

template<typename T>
void maskedSubtract(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                    MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    std::vector<MultiNTensorField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Masked_A_minus_B_NTensor3D<T>, domain, fields, mask);
}

template<typename T>
MultiNTensorField3D<T>* subtract( MultiNTensorField3D<T>& A,
                                  MultiNTensorField3D<T>& B, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<T>* result =
        generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    subtract(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<T>* maskedSubtract( MultiNTensorField3D<T>& A,
                                        MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<T>* result =
        generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    maskedSubtract(A, B, *result, mask, domain);
    return result;
}


template<typename T>
void subtract(MultiNTensorField3D<T>& A, T* alpha, int size,
              MultiNTensorField3D<T>& result, Box3D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_minus_alpha_NTensor3D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedSubtract(MultiNTensorField3D<T>& A, T* alpha, int size,
                    MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_minus_alpha_NTensor3D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField3D<T>* subtract( MultiNTensorField3D<T>& A,
                                  T* alpha, int size, Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    subtract(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<T>* maskedSubtract( MultiNTensorField3D<T>& A,
                                        T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedSubtract(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void subtract(T* alpha, int size, MultiNTensorField3D<T>& A,
              MultiNTensorField3D<T>& result, Box3D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Alpha_minus_A_NTensor3D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedSubtract(T* alpha, int size, MultiNTensorField3D<T>& A,
                    MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_Alpha_minus_A_NTensor3D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField3D<T>* subtract(T* alpha, int size,
                                 MultiNTensorField3D<T>& A, Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    subtract(alpha, size, A, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<T>* maskedSubtract(T* alpha, int size,
                                       MultiNTensorField3D<T>& A, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedSubtract(alpha, size, A, *result, mask, domain);
    return result;
}

template<typename T>
void multiply(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
              MultiNTensorField3D<T>& result, Box3D domain)
{
    std::vector<MultiNTensorField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_times_B_NTensor3D<T>, domain, fields );
}

template<typename T>
void maskedMultiply(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                    MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    std::vector<MultiNTensorField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Masked_A_times_B_NTensor3D<T>, domain, fields, mask );
}

template<typename T>
MultiNTensorField3D<T>* multiply(MultiNTensorField3D<T>& A,
                                 MultiNTensorField3D<T>& B, Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    multiply(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<T>* maskedMultiply(MultiNTensorField3D<T>& A,
                                       MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedMultiply(A, B, *result, mask, domain);
    return result;
}


template<typename T>
void multiply(MultiNTensorField3D<T>& A, T* alpha, int size,
              MultiNTensorField3D<T>& result, Box3D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_times_alpha_NTensor3D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedMultiply(MultiNTensorField3D<T>& A, T* alpha, int size,
                    MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_times_alpha_NTensor3D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField3D<T>* multiply(MultiNTensorField3D<T>& A,
                                 T* alpha, int size, Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    multiply(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<T>* maskedMultiply(MultiNTensorField3D<T>& A,
                                       T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedMultiply(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void divide(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
            MultiNTensorField3D<T>& result, Box3D domain)
{
    std::vector<MultiNTensorField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_dividedBy_B_NTensor3D<T>, domain, fields );
}

template<typename T>
void maskedDivide(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                  MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    std::vector<MultiNTensorField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Masked_A_dividedBy_B_NTensor3D<T>, domain, fields, mask );
}

template<typename T>
MultiNTensorField3D<T>* divide(MultiNTensorField3D<T>& A,
                               MultiNTensorField3D<T>& B, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<T>* result =
        generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    divide(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<T>* maskedDivide(MultiNTensorField3D<T>& A,
                                     MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<T>* result =
        generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    maskedDivide(A, B, *result, mask, domain);
    return result;
}


template<typename T>
void divide(MultiNTensorField3D<T>& A, T* alpha, int size,
            MultiNTensorField3D<T>& result, Box3D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_dividedBy_alpha_NTensor3D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedDivide(MultiNTensorField3D<T>& A, T* alpha, int size,
                  MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_dividedBy_alpha_NTensor3D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField3D<T>* divide(MultiNTensorField3D<T>& A,
                               T* alpha, int size, Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    divide(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<T>* maskedDivide(MultiNTensorField3D<T>& A,
                                     T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedDivide(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void divide(T* alpha, int size, MultiNTensorField3D<T>& A,
            MultiNTensorField3D<T>& result, Box3D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Alpha_dividedBy_A_NTensor3D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedDivide(T* alpha, int size, MultiNTensorField3D<T>& A,
                  MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_Alpha_dividedBy_A_NTensor3D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField3D<T>* divide(T* alpha, int size, MultiNTensorField3D<T>& A,
                               Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    divide(alpha, size, A, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<T>* maskedDivide(T* alpha, int size, MultiNTensorField3D<T>& A, MultiNTensorField3D<int>& mask,
                                     Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedDivide(alpha, size, A, *result, mask, domain);
    return result;
}


template<typename T>
void toThePower(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                MultiNTensorField3D<T>& result, Box3D domain)
{
    std::vector<MultiNTensorField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_toThePower_B_NTensor3D<T>, domain, fields );
}

template<typename T>
void maskedToThePower(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                      MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    std::vector<MultiNTensorField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Masked_A_toThePower_B_NTensor3D<T>, domain, fields, mask );
}

template<typename T>
MultiNTensorField3D<T>* toThePower(MultiNTensorField3D<T>& A,
                                   MultiNTensorField3D<T>& B, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<T>* result =
        generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    toThePower(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<T>* maskedToThePower(MultiNTensorField3D<T>& A,
                                         MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<T>* result =
        generateIntersectMultiNTensorField<T>(A,B, domain, A.getNdim());
    maskedToThePower(A, B, *result, mask, domain);
    return result;
}


template<typename T>
void toThePower(MultiNTensorField3D<T>& A, T* alpha, int size,
                MultiNTensorField3D<T>& result, Box3D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_toThePower_alpha_NTensor3D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedToThePower(MultiNTensorField3D<T>& A, T* alpha, int size,
                      MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_toThePower_alpha_NTensor3D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField3D<T>* toThePower(MultiNTensorField3D<T>& A,
                                   T* alpha, int size, Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    toThePower(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<T>* maskedToThePower(MultiNTensorField3D<T>& A,
                                         T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedToThePower(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void toThePower(T* alpha, int size, MultiNTensorField3D<T>& A,
                MultiNTensorField3D<T>& result, Box3D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Alpha_toThePower_A_NTensor3D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedToThePower(T* alpha, int size, MultiNTensorField3D<T>& A,
                      MultiNTensorField3D<T>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_Alpha_toThePower_A_NTensor3D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField3D<T>* toThePower(T* alpha, int size, MultiNTensorField3D<T>& A,
                                   Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    toThePower(alpha, size, A, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<T>* maskedToThePower(T* alpha, int size, MultiNTensorField3D<T>& A, MultiNTensorField3D<int>& mask,
                                         Box3D domain)
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(A, domain, A.getNdim());
    maskedToThePower(alpha, size, A, *result, mask, domain);
    return result;
}


template<typename T>
void equals(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
            MultiNTensorField3D<int>& result, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_equals_B_NTensor3D<T>, domain, fields );
}

template<typename T>
void maskedEquals(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                  MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_A_equals_B_NTensor3D<T>, domain, fields );
}

template<typename T>
MultiNTensorField3D<int>* equals(MultiNTensorField3D<T>& A,
                                 MultiNTensorField3D<T>& B, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    equals(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<int>* maskedEquals(MultiNTensorField3D<T>& A,
                                       MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    maskedEquals(A, B, *result, mask, domain);
    return result;
}

template<typename T>
void equals( MultiNTensorField3D<T>& A, T* alpha, int size,
             MultiNTensorField3D<int>& result, Box3D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_equals_alpha_NTensor3D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedEquals( MultiNTensorField3D<T>& A, T* alpha, int size,
                   MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_equals_alpha_NTensor3D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField3D<int>* equals(MultiNTensorField3D<T>& A,
                                 T* alpha, int size, Box3D domain)
{
    MultiNTensorField3D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    equals(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<int>* maskedEquals(MultiNTensorField3D<T>& A,
                                       T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    maskedEquals(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void lessThan(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
              MultiNTensorField3D<int>& result, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_lessThan_B_NTensor3D<T>, domain, fields );
}

template<typename T>
void maskedLessThan(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                    MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_A_lessThan_B_NTensor3D<T>, domain, fields );
}

template<typename T>
MultiNTensorField3D<int>* lessThan(MultiNTensorField3D<T>& A,
                                   MultiNTensorField3D<T>& B, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    lessThan(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<int>* maskedLessThan(MultiNTensorField3D<T>& A,
                                         MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    maskedLessThan(A, B, *result, mask, domain);
    return result;
}

template<typename T>
void lessThan( MultiNTensorField3D<T>& A, T* alpha, int size,
               MultiNTensorField3D<int>& result, Box3D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_lessThan_alpha_NTensor3D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedLessThan( MultiNTensorField3D<T>& A, T* alpha, int size,
                     MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_lessThan_alpha_NTensor3D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField3D<int>* lessThan(MultiNTensorField3D<T>& A,
                                   T* alpha, int size, Box3D domain)
{
    MultiNTensorField3D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    lessThan(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<int>* maskedLessThan(MultiNTensorField3D<T>& A,
                                         T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    maskedLessThan(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void lessEqual(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
               MultiNTensorField3D<int>& result, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_lessEqual_B_NTensor3D<T>, domain, fields );
}

template<typename T>
void maskedLessEqual(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                     MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_A_lessEqual_B_NTensor3D<T>, domain, fields );
}

template<typename T>
MultiNTensorField3D<int>* lessEqual(MultiNTensorField3D<T>& A,
                                    MultiNTensorField3D<T>& B, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    lessEqual(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<int>* maskedLessEqual(MultiNTensorField3D<T>& A,
                                          MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    maskedLessEqual(A, B, *result, mask, domain);
    return result;
}

template<typename T>
void lessEqual( MultiNTensorField3D<T>& A, T* alpha, int size,
                MultiNTensorField3D<int>& result, Box3D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_lessEqual_alpha_NTensor3D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedLessEqual( MultiNTensorField3D<T>& A, T* alpha, int size,
                      MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_lessEqual_alpha_NTensor3D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField3D<int>* lessEqual(MultiNTensorField3D<T>& A,
                                    T* alpha, int size, Box3D domain)
{
    MultiNTensorField3D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    lessEqual(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<int>* maskedLessEqual(MultiNTensorField3D<T>& A,
                                          T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    maskedLessEqual(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void greaterThan(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                 MultiNTensorField3D<int>& result, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_greaterThan_B_NTensor3D<T>, domain, fields );
}

template<typename T>
void maskedGreaterThan(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                       MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_A_greaterThan_B_NTensor3D<T>, domain, fields );
}

template<typename T>
MultiNTensorField3D<int>* greaterThan(MultiNTensorField3D<T>& A,
                                      MultiNTensorField3D<T>& B, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    greaterThan(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<int>* maskedGreaterThan(MultiNTensorField3D<T>& A,
                                            MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    maskedGreaterThan(A, B, *result, mask, domain);
    return result;
}

template<typename T>
void greaterThan( MultiNTensorField3D<T>& A, T* alpha, int size,
                  MultiNTensorField3D<int>& result, Box3D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_greaterThan_alpha_NTensor3D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedGreaterThan( MultiNTensorField3D<T>& A, T* alpha, int size,
                        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_greaterThan_alpha_NTensor3D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField3D<int>* greaterThan(MultiNTensorField3D<T>& A,
                                      T* alpha, int size, Box3D domain)
{
    MultiNTensorField3D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    greaterThan(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<int>* maskedGreaterThan(MultiNTensorField3D<T>& A,
                                      T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    maskedGreaterThan(A, alpha, size, *result, mask, domain);
    return result;
}


template<typename T>
void greaterEqual(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                  MultiNTensorField3D<int>& result, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_greaterEqual_B_NTensor3D<T>, domain, fields );
}

template<typename T>
void maskedGreaterEqual(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                        MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_A_greaterEqual_B_NTensor3D<T>, domain, fields );
}

template<typename T>
MultiNTensorField3D<int>* greaterEqual(MultiNTensorField3D<T>& A,
                                       MultiNTensorField3D<T>& B, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    greaterEqual(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<int>* maskedGreaterEqual(MultiNTensorField3D<T>& A,
                                             MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    maskedGreaterEqual(A, B, *result, mask, domain);
    return result;
}

template<typename T>
void greaterEqual( MultiNTensorField3D<T>& A, T* alpha, int size,
                   MultiNTensorField3D<int>& result, Box3D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_greaterEqual_alpha_NTensor3D<T>(alphaVect), domain, A, result );
}

template<typename T>
void maskedGreaterEqual( MultiNTensorField3D<T>& A, T* alpha, int size,
                         MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain )
{
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_greaterEqual_alpha_NTensor3D<T>(alphaVect), domain, A, result, mask );
}

template<typename T>
MultiNTensorField3D<int>* greaterEqual(MultiNTensorField3D<T>& A,
                                       T* alpha, int size, Box3D domain)
{
    MultiNTensorField3D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    greaterEqual(A, alpha, size, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<int>* maskedGreaterEqual(MultiNTensorField3D<T>& A,
                                             T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    maskedGreaterEqual(A, alpha, size, *result, mask, domain);
    return result;
}

template<typename T>
void logicalAnd(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                MultiNTensorField3D<int>& result, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_and_B_NTensor3D<T>, domain, fields );
}

template<typename T>
void maskedLogicalAnd(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                      MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_A_and_B_NTensor3D<T>, domain, fields );
}

template<typename T>
MultiNTensorField3D<int>* logicalAnd(MultiNTensorField3D<T>& A,
                                     MultiNTensorField3D<T>& B, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    logicalAnd(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<int>* maskedLogicalAnd(MultiNTensorField3D<T>& A,
                                           MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    maskedLogicalAnd(A, B, *result, mask, domain);
    return result;
}


template<typename T>
void logicalOr(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
               MultiNTensorField3D<int>& result, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_or_B_NTensor3D<T>, domain, fields );
}

template<typename T>
void maskedLogicalOr(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B,
                     MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_A_or_B_NTensor3D<T>, domain, fields );
}

template<typename T>
MultiNTensorField3D<int>* logicalOr(MultiNTensorField3D<T>& A,
                                    MultiNTensorField3D<T>& B, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    logicalOr(A, B, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<int>* maskedLogicalOr(MultiNTensorField3D<T>& A,
                                          MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    MultiNTensorField3D<int>* result =
        generateIntersectMultiNTensorField<int>(A,B, domain, A.getNdim());
    maskedLogicalOr(A, B, *result, mask, domain);
    return result;
}

template<typename T>
void negate( MultiNTensorField3D<T>& A,
             MultiNTensorField3D<int>& result, Box3D domain )
{
    applyProcessingFunctional (
            new Not_A_NTensor3D<T>(), domain, A, result );
}

template<typename T>
void maskedNegate( MultiNTensorField3D<T>& A,
                   MultiNTensorField3D<int>& result, MultiNTensorField3D<int>& mask, Box3D domain )
{
    applyProcessingFunctional (
            new Masked_Not_A_NTensor3D<T>(), domain, A, result, mask );
}

template<typename T>
MultiNTensorField3D<int>* negate(MultiNTensorField3D<T>& A, Box3D domain)
{
    MultiNTensorField3D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    negate(A, *result, domain);
    return result;
}

template<typename T>
MultiNTensorField3D<int>* maskedNegate(MultiNTensorField3D<T>& A, MultiNTensorField3D<int>& mask, Box3D domain)
{
    MultiNTensorField3D<int>* result =
        generateMultiNTensorField<int>(A, domain, A.getNdim());
    maskedNegate(A, *result, mask, domain);
    return result;
}


/* *************** NTensorField - NTensorField inplace operations *************** */

template<typename T>
void assignInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, Box3D domain) {
    applyProcessingFunctional (
            new AssignNTensorFunctional3D<T>, domain, A, B );
}

template<typename T>
void maskedAssignInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain) {
    applyProcessingFunctional (
            new MaskedAssignNTensorFunctional3D<T>, domain, A, B, mask );
}


template<typename T>
void assignInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, Box3D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new IniConstNTensorFunctional3D<T>(alphaVect), domain, A );
}

template<typename T>
void maskedAssignInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new MaskedIniConstNTensorFunctional3D<T>(alphaVect), domain, A, mask );
}


template<typename T>
void addInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, Box3D domain) {
    applyProcessingFunctional (
            new A_plus_B_inplace_NTensor3D<T>, domain, A, B );
}

template<typename T>
void maskedAddInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain) {
    applyProcessingFunctional (
            new Masked_A_plus_B_inplace_NTensor3D<T>, domain, A, B, mask );
}


template<typename T>
void addInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, Box3D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_plus_alpha_inplace_NTensor3D<T>(alphaVect), domain, A );
}

template<typename T>
void maskedAddInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_plus_alpha_inplace_NTensor3D<T>(alphaVect), domain, A, mask );
}


template<typename T>
void subtractInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, Box3D domain) {
    applyProcessingFunctional (
            new A_minus_B_inplace_NTensor3D<T>, domain, A, B );
}

template<typename T>
void maskedSubtractInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain) {
    applyProcessingFunctional (
            new Masked_A_minus_B_inplace_NTensor3D<T>, domain, A, B, mask );
}


template<typename T>
void subtractInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, Box3D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_minus_alpha_inplace_NTensor3D<T>(alphaVect), domain, A );
}

template<typename T>
void maskedSubtractInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_minus_alpha_inplace_NTensor3D<T>(alphaVect), domain, A, mask );
}


template<typename T>
void multiplyInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, Box3D domain) {
    applyProcessingFunctional (
            new A_times_B_inplace_NTensor3D<T>, domain, A, B );
}

template<typename T>
void maskedMultiplyInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain) {
    applyProcessingFunctional (
            new Masked_A_times_B_inplace_NTensor3D<T>, domain, A, B, mask );
}


template<typename T>
void multiplyInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, Box3D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_times_alpha_inplace_NTensor3D<T>(alphaVect), domain, A );
}

template<typename T>
void maskedMultiplyInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_times_alpha_inplace_NTensor3D<T>(alphaVect), domain, A, mask );
}


template<typename T>
void divideInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, Box3D domain) {
    applyProcessingFunctional (
            new A_dividedBy_B_inplace_NTensor3D<T>, domain, A, B );
}

template<typename T>
void maskedDivideInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain) {
    applyProcessingFunctional (
            new Masked_A_dividedBy_B_inplace_NTensor3D<T>, domain, A, B, mask );
}


template<typename T>
void divideInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, Box3D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_dividedBy_alpha_inplace_NTensor3D<T>(alphaVect), domain, A );
}

template<typename T>
void maskedDivideInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_dividedBy_alpha_inplace_NTensor3D<T>(alphaVect), domain, A, mask );
}


template<typename T>
void toThePowerInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, Box3D domain) {
    applyProcessingFunctional (
            new A_toThePower_B_inplace_NTensor3D<T>, domain, A, B );
}

template<typename T>
void maskedToThePowerInPlace(MultiNTensorField3D<T>& A, MultiNTensorField3D<T>& B, MultiNTensorField3D<int>& mask, Box3D domain) {
    applyProcessingFunctional (
            new Masked_A_toThePower_B_inplace_NTensor3D<T>, domain, A, B, mask );
}


template<typename T>
void toThePowerInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, Box3D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new A_toThePower_alpha_inplace_NTensor3D<T>(alphaVect), domain, A );
}

template<typename T>
void maskedToThePowerInPlace(MultiNTensorField3D<T>& A, T* alpha, int size, MultiNTensorField3D<int>& mask, Box3D domain) {
    PLB_PRECONDITION( size == (int)A.getNdim() );
    std::vector<T> alphaVect(alpha, alpha+size);
    applyProcessingFunctional (
            new Masked_A_toThePower_alpha_inplace_NTensor3D<T>(alphaVect), domain, A, mask );
}

/* *************** UPO ******************************* */

template<typename T>
T compute_UPO_ScalarProduct(MultiNTensorField3D<T>& a, MultiNTensorField3D<T>& b,
                            Box3D domain)
{
    PLB_PRECONDITION(a.getNdim() == b.getNdim());
    UPO_ScalarProductFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, a,b);
    return functional.getSum();
}

template<typename T>
T masked_compute_UPO_ScalarProduct(MultiNTensorField3D<T>& a, MultiNTensorField3D<T>& b,
                                   MultiNTensorField3D<int>& mask, Box3D domain)
{
    PLB_PRECONDITION(a.getNdim() == b.getNdim());
    Masked_UPO_ScalarProductFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, a,b,mask);
    return functional.getSum();
}

}  // namespace plb

#endif  // SWIG_DATA_ANALYSIS_WRAPPER_3D_HH
