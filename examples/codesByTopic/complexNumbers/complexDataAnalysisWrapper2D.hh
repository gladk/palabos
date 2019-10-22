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
#ifndef COMPLEX_DATA_ANALYSIS_WRAPPER_2D_HH
#define COMPLEX_DATA_ANALYSIS_WRAPPER_2D_HH

#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "typeConverterFunctional2D.hh"


namespace plb {

/* ******************************************************************* */
/* *************** PART I. Atomic-block wrappers: Block-Lattice ****** */
/* ******************************************************************* */

/* *************** From Complex to real MultiScalar and MultiTensor Fields ******************************* */

template<typename T, typename U>
void realPart( ScalarField2D<T>& complexField,
               ScalarField2D<U>& realField,
               Box2D domain) {
    applyProcessingFunctional(new FromComplexToRealScalarFieldFunctional2D<T,U> , domain, complexField, realField);
}

template<typename T, typename U>
std::auto_ptr<ScalarField2D<U> > realPart(ScalarField2D<T>& complexField, Box2D domain) {
    ScalarField2D<U> *realField = new ScalarField2D<U>(complexField.getNx(),complexField.getNy());
    realPart(complexField, *realField, domain);
    return std::auto_ptr<ScalarField2D<U> >(realField);
}

template<typename T, typename U>
void imaginaryPart( ScalarField2D<T>& complexField,
                    ScalarField2D<U>& realField,
                    Box2D domain) {
    applyProcessingFunctional(new FromComplexToImaginaryScalarFieldFunctional2D<T,U> , domain, complexField, realField);
}

template<typename T, typename U>
std::auto_ptr<ScalarField2D<U> > imaginaryPart(ScalarField2D<T>& complexField, Box2D domain) {
    ScalarField2D<U> *realField = new ScalarField2D<U>(complexField.getNx(),complexField.getNy());
    imaginaryPart(complexField, *realField, domain);
    return std::auto_ptr<ScalarField2D<U> >(realField);
}

/* ***** For TensorFields ****/    


template<typename T, typename R, int d>
void realPart( TensorField2D<T,d>& complexField,
               TensorField2D<R,d>& realField,
               Box2D domain) {
    applyProcessingFunctional(new FromComplexToRealTensorFieldFunctional2D<T,R,d> , domain, complexField, realField);
}

template<typename T, typename R, int d>
std::auto_ptr<TensorField2D<R,d> > realPart(TensorField2D<T,d>& complexField, Box2D domain) {
    TensorField2D<R,d> *realField = new TensorField2D<R,d>(complexField.getNx(),complexField.getNy());
    realPart(complexField, *realField, domain);
    return std::auto_ptr<TensorField2D<R,d> >(realField);
}


template<typename T, typename R, int d>
void imaginaryPart(TensorField2D<T,d>& complexField,
                   TensorField2D<R,d>& realField,
                   Box2D domain) {
    applyProcessingFunctional(new FromComplexToImaginaryTensorFieldFunctional2D<T,R,d> , domain, complexField, realField);
}

template<typename T, typename R, int d>
std::auto_ptr<TensorField2D<R,d> > imaginaryPart(TensorField2D<T,d>& complexField, Box2D domain) {
    TensorField2D<R,d> *realField = new TensorField2D<R,d>(complexField.getNx(),complexField.getNy());
    imaginaryPart(complexField, *realField, domain);
    return std::auto_ptr<TensorField2D<R,d> >(realField);
}   

/* ************************************************************************* */
/* *************** PART IV. Multi-block wrappers: Multi-Block-Lattice ****** */
/* ************************************************************************* */

/* *************** From Complex to real MultiScalar and MultiTensor Fields ******************************* */

template<typename T, typename U>
void realPart( MultiScalarField2D<T>& complexField,
               MultiScalarField2D<U>& realField,
               Box2D domain) {
    applyProcessingFunctional(new FromComplexToRealScalarFieldFunctional2D<T,U> , domain, complexField, realField);
}

template<typename T, typename U>
std::auto_ptr<MultiScalarField2D<U> > realPart(MultiScalarField2D<T>& complexField, Box2D domain) {
    MultiScalarField2D<U> *realField = new MultiScalarField2D<U>(complexField.getNx(),complexField.getNy());
    realPart(complexField, *realField, domain);
    return std::auto_ptr<MultiScalarField2D<U> >(realField);
}

template<typename T, typename U>
void imaginaryPart( MultiScalarField2D<T>& complexField,
                    MultiScalarField2D<U>& realField,
                    Box2D domain) {
    applyProcessingFunctional(new FromComplexToImaginaryScalarFieldFunctional2D<T,U> , domain, complexField, realField);
}

template<typename T, typename U>
std::auto_ptr<MultiScalarField2D<U> > imaginaryPart(MultiScalarField2D<T>& complexField, Box2D domain) {
    MultiScalarField2D<U> *realField = new MultiScalarField2D<U>(complexField.getNx(),complexField.getNy());
    imaginaryPart(complexField, *realField, domain);
    return std::auto_ptr<MultiScalarField2D<U> >(realField);
}

/* ***** For MultiTensorFields ****/    


template<typename T, typename R, int d>
void realPart( MultiTensorField2D<T,d>& complexField,
               MultiTensorField2D<R,d>& realField,
               Box2D domain) {
    applyProcessingFunctional(new FromComplexToRealTensorFieldFunctional2D<T,R,d> , domain, complexField, realField);
}

template<typename T, typename R, int d>
std::auto_ptr<MultiTensorField2D<R,d> > realPart(MultiTensorField2D<T,d>& complexField, Box2D domain) {
    MultiTensorField2D<R,d> *realField = new MultiTensorField2D<R,d>(complexField.getNx(),complexField.getNy());
    realPart(complexField, *realField, domain);
    return std::auto_ptr<MultiTensorField2D<R,d> >(realField);
}


template<typename T, typename R, int d>
void imaginaryPart(MultiTensorField2D<T,d>& complexField,
                   MultiTensorField2D<R,d>& realField,
                   Box2D domain) {
    applyProcessingFunctional(new FromComplexToImaginaryTensorFieldFunctional2D<T,R,d> , domain, complexField, realField);
}

template<typename T, typename R, int d>
std::auto_ptr<MultiTensorField2D<R,d> > imaginaryPart(MultiTensorField2D<T,d>& complexField, Box2D domain) {
    MultiTensorField2D<R,d> *realField = new MultiTensorField2D<R,d>(complexField.getNx(),complexField.getNy());
    imaginaryPart(complexField, *realField, domain);
    return std::auto_ptr<MultiTensorField2D<R,d> >(realField);
}   

}  // namespace plb

#endif  // DATA_ANALYSIS_WRAPPER_2D_HH
