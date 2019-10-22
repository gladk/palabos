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
 * Helper functions for data analysis -- header file.
 */

#ifndef ML_DATA_ANALYSIS_WRAPPER_3D_HH
#define ML_DATA_ANALYSIS_WRAPPER_3D_HH

#include "dataProcessors/dataAnalysisFunctional3D.h"
#include "gridRefinement/dataAnalysisWrapper3D.h"
#include "gridRefinement/multiLevelWrapper3D.h"
#include "gridRefinement/multiLevelFieldGenerator3D.h"

namespace plb {

/* ******************************************************************* */
/* *************** PART I. Multi-block wrappers: Block-Lattice ****** */
/* ******************************************************************* */

/* *************** Density ****************************************** */

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
void computeDensity(
    MultiLevelActions3D<T,Descriptor,Engine>& lattices, 
    MultiLevelScalarField3D<T>& densities, 
    Box3D domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new BoxDensityFunctional3D<T,Descriptor>(), domain, levelOfDomain,
            lattices, densities );
}

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelScalarField3D<T> > computeDensity(
    MultiLevelActions3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelScalarField3D<T> > densities 
        = generateMultiLevelScalarField3D<T>(lattices.getOgs(), domain, levelOfDomain);

    computeDensity(lattices, *densities, domain, levelOfDomain);
    
    return densities;
}

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelScalarFieldForOutput3D<T> > computeDensity(
    MultiLevelActions3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain, bool crop)
{
    return exportForOutput(*computeDensity(lattices, domain, levelOfDomain), domain, levelOfDomain, crop);
}

/* *************** Velocity ****************************************** */

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
void computeVelocity(
    MultiLevelActions3D<T,Descriptor,Engine>& lattices, 
    MultiLevelTensorField3D<T,Descriptor<T>::d> &velocities, 
    Box3D domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new BoxVelocityFunctional3D<T,Descriptor>(), domain, levelOfDomain,
            lattices, velocities );
}

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelTensorField3D<T,Descriptor<T>::d> > computeVelocity(
    MultiLevelActions3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelTensorField3D<T,Descriptor<T>::d> > velocities = 
        generateMultiLevelTensorField3D<T,Descriptor<T>::d>(lattices.getOgs(), domain, levelOfDomain);

    computeVelocity(lattices, *velocities, domain, levelOfDomain);

    return velocities;
}

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelTensorFieldForOutput3D<T,Descriptor<T>::d> > computeVelocity(
    MultiLevelActions3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain, bool crop)
{
    return exportForOutput(*computeVelocity(lattices, domain, levelOfDomain), domain, levelOfDomain, crop);
}





/* ******************************************************************* */
/* *************** PART I. Multi-block wrappers: Block-Lattice ****** */
/* ******************************************************************* */

/* *************** Density ****************************************** */

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
void computeDensity(
    MultiLevelCoupling3D<T,Descriptor,Engine>& lattices, 
    MultiLevelScalarField3D<T>& densities, 
    Box3D domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new BoxDensityFunctional3D<T,Descriptor>(), domain, levelOfDomain,
            lattices, densities );
}

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelScalarField3D<T> > computeDensity(
    MultiLevelCoupling3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelScalarField3D<T> > densities 
        = generateMultiLevelScalarField3D<T>(lattices.getOgs(), domain, levelOfDomain);

    computeDensity(lattices, *densities, domain, levelOfDomain);
    
    return densities;
}

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelScalarFieldForOutput3D<T> > computeDensity(
    MultiLevelCoupling3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain, bool crop)
{
    return exportForOutput(*computeDensity(lattices, domain, levelOfDomain), domain, levelOfDomain, crop);
}

/* *************** Velocity ****************************************** */

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
void computeVelocity(
    MultiLevelCoupling3D<T,Descriptor,Engine>& lattices, 
    MultiLevelTensorField3D<T,Descriptor<T>::d> &velocities, 
    Box3D domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new BoxVelocityFunctional3D<T,Descriptor>(), domain, levelOfDomain,
            lattices, velocities );
}

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelTensorField3D<T,Descriptor<T>::d> > computeVelocity(
    MultiLevelCoupling3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelTensorField3D<T,Descriptor<T>::d> > velocities = 
        generateMultiLevelTensorField3D<T,Descriptor<T>::d>(lattices.getOgs(), domain, levelOfDomain);

    computeVelocity(lattices, *velocities, domain, levelOfDomain);

    return velocities;
}

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelTensorFieldForOutput3D<T,Descriptor<T>::d> > computeVelocity(
    MultiLevelCoupling3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain, bool crop)
{
    return exportForOutput(*computeVelocity(lattices, domain, levelOfDomain), domain, levelOfDomain, crop);
}




/* ******************************************************************* */
/* *************** PART II. Multi-block wrappers: Scalar-Field ******** */
/* ******************************************************************* */

/* ========================= copy ============================ */
template<typename T>
void copy(MultiLevelScalarField3D<T> &from, 
          const Box3D &fromDomain, plint levelOfFromDomain,
          MultiLevelScalarFieldForOutput3D<T> &to, 
          const Box3D &toDomain, plint levelOfToDomain)
{
    PLB_ASSERT(from.getNumLevels() == to.getNumLevels() && 
        "MultiScalarTensorField3D (from) must have the save number of levels as the MultiLevelScalarFieldForOutput3D (to)");

    for (plint iA = 0; iA < from.getNumLevels(); ++iA) {
        Box3D currentFromDomain = (iA - levelOfFromDomain >= 0) ? 
            fromDomain.multiply(util::intTwoToThePower(iA - levelOfFromDomain)) : 
            fromDomain.divide(util::intTwoToThePower(levelOfFromDomain-iA));

        Box3D currentToDomain = (iA - levelOfToDomain >= 0) ? 
            toDomain.multiply(util::intTwoToThePower(iA - levelOfToDomain)) : 
            toDomain.divide(util::intTwoToThePower(levelOfToDomain-iA));

        copy(from.getLevel(iA), currentFromDomain, to.getLevel(iA), currentToDomain);
    }   

}

/* *************** exportForOutput ****************************************** */

template<typename T>
std::auto_ptr<MultiLevelScalarFieldForOutput3D<T> > 
  exportForOutput(MultiLevelScalarField3D<T> &from, 
          const Box3D &domain, plint levelOfDomain, bool crop)
{
    std::auto_ptr<MultiLevelScalarFieldForOutput3D<T> > output = 
        generateMultiLevelScalarFieldForOutput3D<T>(from.getOgs(), domain, levelOfDomain, crop);

    copy(from, domain, levelOfDomain, *output, domain, levelOfDomain);

    return output;
}

/* *************** Extract Sub-MultiLevelTensorField *************************** */

template<typename T>
void extractSubDomain( MultiLevelScalarField3D<T> &field,
                       MultiLevelScalarField3D<T> &extractedField,
                       const Box3D &domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new ExtractScalarSubDomainFunctional3D<T>, domain, levelOfDomain, field, extractedField );
}

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> > extractSubDomain(
    MultiLevelScalarField3D<T>& field, 
    const Box3D &domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelScalarField3D<T> > extractedField = 
        generateMultiLevelScalarField3D<T>(field.getOgs(), domain, levelOfDomain);

    extractSubDomain(field, *extractedField, domain, levelOfDomain);

    return extractedField;
}

/* *************** add ****************************************** */

template<typename T>
void add(MultiLevelScalarField3D<T> &A, 
         MultiLevelScalarField3D<T> &B, 
         MultiLevelScalarField3D<T> &aPlusB,
         const Box3D &domain, plint levelOfDomain) 
{
    std::vector<MultiLevelScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&aPlusB);
    applyProcessingFunctional (
            new A_plus_B_functional3D<T>, domain, levelOfDomain, fields );
}

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> > 
    add(MultiLevelScalarField3D<T> &A, 
        MultiLevelScalarField3D<T> &B, 
        const Box3D &domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelScalarField3D<T> > aPlusB = 
        generateMultiLevelScalarField3D<T>(A.getOgs(), domain, levelOfDomain);

    add(A, B, *aPlusB, domain, levelOfDomain);

    return aPlusB;
}

/* *************** addInPlace ****************************************** */

template<typename T>
void addInPlace(MultiLevelScalarField3D<T> &A, 
                MultiLevelScalarField3D<T> &B, 
                const Box3D &domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new A_plus_B_inplace_functional3D<T>, 
            domain, levelOfDomain, A, B );
}

/* *************** subtract ****************************************** */

template<typename T>
void subtract(MultiLevelScalarField3D<T> &A, 
              MultiLevelScalarField3D<T> &B, 
              MultiLevelScalarField3D<T> &aMinusB,
              const Box3D &domain, plint levelOfDomain)
{
    std::vector<MultiLevelScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&aMinusB);
    applyProcessingFunctional (
            new A_minus_B_functional3D<T>, domain, levelOfDomain, fields );
}

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> >  
    subtract(MultiLevelScalarField3D<T> &A, 
             MultiLevelScalarField3D<T> &B, 
             const Box3D &domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelScalarField3D<T> > aMinusB = 
        generateMultiLevelScalarField3D<T>(A.getOgs(), domain, levelOfDomain);

    subtract(A, B, *aMinusB, domain, levelOfDomain);

    return aMinusB;
}

/* *************** subtractInPlace ****************************************** */

template<typename T>
void subtractInPlace(MultiLevelScalarField3D<T> &A, 
              MultiLevelScalarField3D<T> &B, 
              const Box3D &domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new A_minus_B_inplace_functional3D<T>, 
            domain, levelOfDomain, A, B );
}

/* *************** multiply ****************************************** */

template<typename T>
void multiply(MultiLevelScalarField3D<T> &A, 
              MultiLevelScalarField3D<T> &B, 
              MultiLevelScalarField3D<T> &aTimesB,
              const Box3D &domain, plint levelOfDomain)
{
    std::vector<MultiLevelScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&aTimesB);
    applyProcessingFunctional (
            new A_times_B_functional3D<T>, domain, levelOfDomain, fields );
}

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> >   
    multiply(MultiLevelScalarField3D<T> &A, 
             MultiLevelScalarField3D<T> &B, 
             const Box3D &domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelScalarField3D<T> > aTimesB = 
        generateMultiLevelScalarField3D<T>(A.getOgs(), domain, levelOfDomain);

    multiply(A, B, *aTimesB, domain, levelOfDomain);

    return aTimesB;
}

/* *************** multiplyInPlace ****************************************** */

template<typename T>
void multiplyInPlace(MultiLevelScalarField3D<T> &A, 
                     MultiLevelScalarField3D<T> &B, 
                     const Box3D &domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new A_times_B_inplace_functional3D<T>, 
            domain, levelOfDomain, A, B );
}

/* *************** multiply by scalar ****************************************** */

template<typename T>
void multiply(MultiLevelScalarField3D<T> &A, 
              T alpha,
              MultiLevelScalarField3D<T> &aTimesAlpha,
              const Box3D &domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new A_times_alpha_functional3D<T>(alpha), domain, levelOfDomain, A, aTimesAlpha );
}

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> >   
    multiply(MultiLevelScalarField3D<T> &A, 
             T alpha,
             const Box3D &domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelScalarField3D<T> > aTimesAlpha = 
        generateMultiLevelScalarField3D<T>(A.getOgs(), domain, levelOfDomain);

    multiply(A, alpha, *aTimesAlpha, domain, levelOfDomain);

    return aTimesAlpha;
}

/* *************** multiplyInPlace by scalar ****************************************** */

template<typename T>
void multiplyInPlace(MultiLevelScalarField3D<T> &A, 
                     T alpha,
                     const Box3D &domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new A_times_alpha_inplace_functional3D<T>(alpha), 
            domain, levelOfDomain, A );
}

/* *************** divide ****************************************** */

template<typename T>
void divide(MultiLevelScalarField3D<T> &A, 
            MultiLevelScalarField3D<T> &B, 
            MultiLevelScalarField3D<T> &aDivideB,
            const Box3D &domain, plint levelOfDomain)
{
    std::vector<MultiLevelScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&aDivideB);
    applyProcessingFunctional (
            new A_dividedBy_B_functional3D<T>, domain, levelOfDomain, fields );
}

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> >  
    divide(MultiLevelScalarField3D<T> &A, 
           MultiLevelScalarField3D<T> &B, 
           const Box3D &domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelScalarField3D<T> > aDivideB = 
        generateMultiLevelScalarField3D<T>(A.getOgs(), domain, levelOfDomain);

    divide(A, B, *aDivideB, domain, levelOfDomain);

    return aDivideB;
}

/* *************** divide ****************************************** */

template<typename T>
void divideInPlace(MultiLevelScalarField3D<T> &A, 
                   MultiLevelScalarField3D<T> &B, 
                   const Box3D &domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new A_dividedBy_B_inplace_functional3D<T>, 
            domain, levelOfDomain, A, B );
}


/* ******************************************************************* */
/* *************** PART III. Multi-block wrappers: Tensor-field ******* */
/* ******************************************************************* */

/* ========================= copy ============================ */
template<typename T, int nDim>
void copy(MultiLevelTensorField3D<T,nDim> &from, 
          const Box3D &fromDomain, plint levelOfFromDomain,
          MultiLevelTensorFieldForOutput3D<T,nDim> &to, 
          const Box3D &toDomain, plint levelOfToDomain)
{
    PLB_ASSERT(from.getNumLevels() == to.getNumLevels() && 
        "MultiLevelTensorField3D (from) must have the save number of levels as the MultiLevelTensorFieldForOutput3D (to)");

    for (plint iA = 0; iA < from.getNumLevels(); ++iA) {
        Box3D currentFromDomain = (iA - levelOfFromDomain >= 0) ? 
            fromDomain.multiply(util::intTwoToThePower(iA - levelOfFromDomain)) : 
            fromDomain.divide(util::intTwoToThePower(levelOfFromDomain-iA));

        Box3D currentToDomain = (iA - levelOfToDomain >= 0) ? 
            toDomain.multiply(util::intTwoToThePower(iA - levelOfToDomain)) : 
            toDomain.divide(util::intTwoToThePower(levelOfToDomain-iA));

        copy(from.getLevel(iA), currentFromDomain, to.getLevel(iA), currentToDomain);
    }   

}

/* *************** exportForOutput ****************************************** */

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorFieldForOutput3D<T,nDim> > 
  exportForOutput(MultiLevelTensorField3D<T,nDim> &from, 
          const Box3D &domain, plint levelOfDomain, bool crop)
{
    std::auto_ptr<MultiLevelTensorFieldForOutput3D<T,nDim> > output = 
        generateMultiLevelTensorFieldForOutput3D<T,nDim>(from.getOgs(), domain, levelOfDomain, crop);

    copy(from, domain, levelOfDomain, *output, domain, levelOfDomain);

    return output;
}

/* *************** Extract Sub-MultiLevelTensorField *************************** */

template<typename T, int nDim>
void extractSubDomain( MultiLevelTensorField3D<T,nDim> &field,
                       MultiLevelTensorField3D<T,nDim> &extractedField,
                       const Box3D &domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new ExtractTensorSubDomainFunctional3D<T,nDim>, domain, levelOfDomain, field, extractedField );
}

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> > extractSubDomain(
    MultiLevelTensorField3D<T,nDim>& field, 
    const Box3D &domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelTensorField3D<T,nDim> > extractedField = 
        generateMultiLevelTensorField3D<T,nDim>(field.getOgs(), domain, levelOfDomain);

    extractSubDomain(field, *extractedField, domain, levelOfDomain);

    return extractedField;
}

/* ***** Component (multilevel scalar-field) out of a multilevel-tensor-field *** */

template<typename T, int nDim>
void extractComponent(
  MultiLevelTensorField3D<T,nDim>& tensorField, 
  MultiLevelScalarField3D<T>& component, 
  const Box3D &domain, plint levelOfDomain, int iComponent)
{
    applyProcessingFunctional (
            new ExtractTensorComponentFunctional3D<T,nDim>(iComponent), 
            domain, levelOfDomain, component, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<MultiLevelScalarField3D<T> > extractComponent(
  MultiLevelTensorField3D<T,nDim>& tensorField, 
  const Box3D &domain, plint levelOfDomain, int iComponent)
{
    std::auto_ptr<MultiLevelScalarField3D<T> > component =
        generateMultiLevelScalarField3D<T>(tensorField.getOgs(), domain, levelOfDomain);
    extractComponent(tensorField, *component, domain, levelOfDomain, iComponent);
    return component;
}

/* *************** add ****************************************** */

template<typename T, int nDim>
void add(MultiLevelTensorField3D<T,nDim> &A, 
         MultiLevelTensorField3D<T,nDim> &B, 
         MultiLevelTensorField3D<T,nDim> &aPlusB,
         const Box3D &domain, plint levelOfDomain) 
{
    std::vector<MultiLevelTensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&aPlusB);
    applyProcessingFunctional (
            new Tensor_A_plus_B_functional3D<T,nDim>, domain, levelOfDomain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> > 
    add(MultiLevelTensorField3D<T,nDim> &A, 
        MultiLevelTensorField3D<T,nDim> &B, 
        const Box3D &domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelTensorField3D<T,nDim> > aPlusB = 
        generateMultiLevelTensorField3D<T,nDim>(A.getOgs(), domain, levelOfDomain);

    add(A, B, *aPlusB, domain, levelOfDomain);

    return aPlusB;
}

/* *************** addInPlace ****************************************** */

template<typename T, int nDim>
void addInPlace(MultiLevelTensorField3D<T,nDim> &A, 
                MultiLevelTensorField3D<T,nDim> &B, 
                const Box3D &domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new Tensor_A_plus_B_inplace_functional3D<T,nDim>, 
            domain, levelOfDomain, A, B );
}

/* *************** subtract ****************************************** */

template<typename T, int nDim>
void subtract(MultiLevelTensorField3D<T,nDim> &A, 
              MultiLevelTensorField3D<T,nDim> &B, 
              MultiLevelTensorField3D<T,nDim> &aMinusB,
              const Box3D &domain, plint levelOfDomain)
{
    std::vector<MultiLevelTensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&aMinusB);
    applyProcessingFunctional (
            new Tensor_A_minus_B_functional3D<T,nDim>, domain, levelOfDomain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> >  
    subtract(MultiLevelTensorField3D<T,nDim> &A, 
             MultiLevelTensorField3D<T,nDim> &B, 
             const Box3D &domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelTensorField3D<T,nDim> > aMinusB = 
        generateMultiLevelTensorField3D<T,nDim>(A.getOgs(), domain, levelOfDomain);

    subtract(A, B, *aMinusB, domain, levelOfDomain);

    return aMinusB;
}

/* *************** subtractInPlace ****************************************** */

template<typename T, int nDim>
void subtractInPlace(MultiLevelTensorField3D<T,nDim> &A, 
              MultiLevelTensorField3D<T,nDim> &B, 
              const Box3D &domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new Tensor_A_minus_B_inplace_functional3D<T,nDim>, 
            domain, levelOfDomain, A, B );
}

/* *************** multiply ****************************************** */

template<typename T, int nDim>
void multiply(MultiLevelTensorField3D<T,nDim> &A, 
              MultiLevelTensorField3D<T,nDim> &B, 
              MultiLevelTensorField3D<T,nDim> &aTimesB,
              const Box3D &domain, plint levelOfDomain)
{
    std::vector<MultiLevelTensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&aTimesB);
    applyProcessingFunctional (
            new Tensor_A_times_B_functional3D<T,nDim>, domain, levelOfDomain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> >   
    multiply(MultiLevelTensorField3D<T,nDim> &A, 
             MultiLevelTensorField3D<T,nDim> &B, 
             const Box3D &domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelTensorField3D<T,nDim> > aTimesB = 
        generateMultiLevelTensorField3D<T,nDim>(A.getOgs(), domain, levelOfDomain);

    multiply(A, B, *aTimesB, domain, levelOfDomain);

    return aTimesB;
}

/* *************** multiplyInPlace ****************************************** */

template<typename T, int nDim>
void multiplyInPlace(MultiLevelTensorField3D<T,nDim> &A, 
                     MultiLevelTensorField3D<T,nDim> &B, 
                     const Box3D &domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new Tensor_A_times_B_inplace_functional3D<T,nDim>, 
            domain, levelOfDomain, A, B );
}

/* *************** multiply by scalar ****************************************** */

template<typename T, int nDim>
void multiply(MultiLevelTensorField3D<T,nDim> &A, 
              T alpha,
              MultiLevelTensorField3D<T,nDim> &aTimesAlpha,
              const Box3D &domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new Tensor_A_times_alpha_functional3D<T,nDim>(alpha), domain, levelOfDomain, A, aTimesAlpha );
}

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> >   
    multiply(MultiLevelTensorField3D<T,nDim> &A, 
             T alpha,
             const Box3D &domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelTensorField3D<T,nDim> > aTimesAlpha = 
        generateMultiLevelTensorField3D<T,nDim>(A.getOgs(), domain, levelOfDomain);

    multiply(A, alpha, *aTimesAlpha, domain, levelOfDomain);

    return aTimesAlpha;
}

/* *************** multiplyInPlace by scalar ****************************************** */

template<typename T, int nDim>
void multiplyInPlace(MultiLevelTensorField3D<T,nDim> &A, 
                     T alpha,
                     const Box3D &domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new Tensor_A_times_alpha_inplace_functional3D<T,nDim>(alpha), 
            domain, levelOfDomain, A );
}

/* *************** divide ****************************************** */

template<typename T, int nDim>
void divide(MultiLevelTensorField3D<T,nDim> &A, 
            MultiLevelTensorField3D<T,nDim> &B, 
            MultiLevelTensorField3D<T,nDim> &aDivideB,
            const Box3D &domain, plint levelOfDomain)
{
    std::vector<MultiLevelTensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&aDivideB);
    applyProcessingFunctional (
            new Tensor_A_dividedBy_B_functional3D<T,nDim>, domain, levelOfDomain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> >  
    divide(MultiLevelTensorField3D<T,nDim> &A, 
           MultiLevelTensorField3D<T,nDim> &B, 
           const Box3D &domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelTensorField3D<T,nDim> > aDivideB = 
        generateMultiLevelTensorField3D<T,nDim>(A.getOgs(), domain, levelOfDomain);

    divide(A, B, *aDivideB, domain, levelOfDomain);

    return aDivideB;
}

/* *************** divide ****************************************** */

template<typename T, int nDim>
void divideInPlace(MultiLevelTensorField3D<T,nDim> &A, 
                   MultiLevelTensorField3D<T,nDim> &B, 
                   const Box3D &domain, plint levelOfDomain)
{
    applyProcessingFunctional (
            new Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>, 
            domain, levelOfDomain, A, B );
}


/* *************** Vorticity ****************************************** */

template<typename T>
std::auto_ptr<MultiLevelTensorField3D<T,3 > > computeVorticity(
    MultiLevelTensorField3D<T,3>& velocities, Box3D domain, plint levelOfDomain)
{
    std::auto_ptr<MultiLevelTensorField3D<T,3> > vorticity = 
        generateMultiLevelTensorField3D<T,3>(velocities.getOgs(), domain, levelOfDomain);

    applyProcessingFunctional (
            new BoxVorticityFunctional3D<T,3>(), domain, levelOfDomain,
            velocities, *vorticity, 1 );
    
    return vorticity;
}

template<typename T>
std::auto_ptr<MultiLevelTensorFieldForOutput3D<T,3 > > computeVorticity(
    MultiLevelTensorField3D<T,3>& velocities, Box3D domain, plint levelOfDomain, bool crop)
{
    std::auto_ptr<MultiLevelTensorField3D<T,3> > vorticity = 
        computeVorticity(velocities, domain, levelOfDomain);

    std::auto_ptr<MultiLevelTensorFieldForOutput3D<T,3> > outputVorticity
        = generateMultiLevelTensorFieldForOutput3D<T,3>(velocities.getOgs(), domain, levelOfDomain, crop);

    copy<T,3>(*vorticity, domain, levelOfDomain, *outputVorticity, domain, levelOfDomain);
    
    return outputVorticity;
}

}  // namespace plb

#endif  // DATA_ANALYSIS_WRAPPER_3D_HH

