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

#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.hh"

namespace plb {

/* *************** Boxed Data Processor functionals ****************** */

template
void applyProcessingFunctional<PRECOMP_T> (
        BoxProcessingFunctional2D_N<PRECOMP_T>* functional,
        Box2D domain, MultiNTensorField2D<PRECOMP_T>& field );
template
void integrateProcessingFunctional<PRECOMP_T> (
        BoxProcessingFunctional2D_N<PRECOMP_T>* functional,
        Box2D domain, MultiNTensorField2D<PRECOMP_T>& field, plint level );


template
void applyProcessingFunctional<PRECOMP_T> (
        BoxProcessingFunctional2D_S<PRECOMP_T>* functional,
        Box2D domain, MultiScalarField2D<PRECOMP_T>& field );
template
void integrateProcessingFunctional<PRECOMP_T> (
        BoxProcessingFunctional2D_S<PRECOMP_T>* functional,
        Box2D domain, MultiScalarField2D<PRECOMP_T>& field, plint level );


template
void applyProcessingFunctional<PRECOMP_T> (
        MaskedBoxProcessingFunctional2D_N<PRECOMP_T>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field,
        MultiNTensorField2D<int>& mask );
template
void integrateProcessingFunctional<PRECOMP_T> (
        MaskedBoxProcessingFunctional2D_N<PRECOMP_T>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field,
        MultiNTensorField2D<int>& mask,
        plint level );

template
void applyProcessingFunctional<PRECOMP_T,int> (
        BoxProcessingFunctional2D_NN<PRECOMP_T,int>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field1,
        MultiNTensorField2D<int>& field2 );
template
void integrateProcessingFunctional<PRECOMP_T,int> (
        BoxProcessingFunctional2D_NN<PRECOMP_T,int>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field1,
        MultiNTensorField2D<int>& field2, plint level );

template
void applyProcessingFunctional<PRECOMP_T,int> (
        BoxProcessingFunctional2D_SS<PRECOMP_T,int>* functional,
        Box2D domain,
        MultiScalarField2D<PRECOMP_T>& field1,
        MultiScalarField2D<int>& field2 );
template
void integrateProcessingFunctional<PRECOMP_T,int> (
        BoxProcessingFunctional2D_SS<PRECOMP_T,int>* functional,
        Box2D domain,
        MultiScalarField2D<PRECOMP_T>& field1,
        MultiScalarField2D<int>& field2, plint level );

template
void applyProcessingFunctional<PRECOMP_T,float> (
        BoxProcessingFunctional2D_NN<PRECOMP_T,float>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field1,
        MultiNTensorField2D<float>& field2 );
template
void integrateProcessingFunctional<PRECOMP_T,float> (
        BoxProcessingFunctional2D_NN<PRECOMP_T,float>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field1,
        MultiNTensorField2D<float>& field2, plint level );

template
void applyProcessingFunctional<PRECOMP_T,float> (
        BoxProcessingFunctional2D_SS<PRECOMP_T,float>* functional,
        Box2D domain,
        MultiScalarField2D<PRECOMP_T>& field1,
        MultiScalarField2D<float>& field2 );
template
void integrateProcessingFunctional<PRECOMP_T,float> (
        BoxProcessingFunctional2D_SS<PRECOMP_T,float>* functional,
        Box2D domain,
        MultiScalarField2D<PRECOMP_T>& field1,
        MultiScalarField2D<float>& field2, plint level );

template
void applyProcessingFunctional<PRECOMP_T,double> (
        BoxProcessingFunctional2D_NN<PRECOMP_T,double>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field1,
        MultiNTensorField2D<double>& field2 );
template
void integrateProcessingFunctional<PRECOMP_T,double> (
        BoxProcessingFunctional2D_NN<PRECOMP_T,double>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field1,
        MultiNTensorField2D<double>& field2, plint level );

template
void applyProcessingFunctional<PRECOMP_T,double> (
        BoxProcessingFunctional2D_SS<PRECOMP_T,double>* functional,
        Box2D domain,
        MultiScalarField2D<PRECOMP_T>& field1,
        MultiScalarField2D<double>& field2 );
template
void integrateProcessingFunctional<PRECOMP_T,double> (
        BoxProcessingFunctional2D_SS<PRECOMP_T,double>* functional,
        Box2D domain,
        MultiScalarField2D<PRECOMP_T>& field1,
        MultiScalarField2D<double>& field2, plint level );

template
void applyProcessingFunctional<PRECOMP_T,float> (
        MaskedBoxProcessingFunctional2D_NN<PRECOMP_T,float>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field1,
        MultiNTensorField2D<float>& field2,
        MultiNTensorField2D<int>& mask );
template
void applyProcessingFunctional<PRECOMP_T,double> (
        MaskedBoxProcessingFunctional2D_NN<PRECOMP_T,double>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field1,
        MultiNTensorField2D<double>& field2,
        MultiNTensorField2D<int>& mask );
template
void applyProcessingFunctional<PRECOMP_T,int> (
        MaskedBoxProcessingFunctional2D_NN<PRECOMP_T,int>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field1,
        MultiNTensorField2D<int>& field2,
        MultiNTensorField2D<int>& mask );
template
void integrateProcessingFunctional<PRECOMP_T,PRECOMP_T> (
        MaskedBoxProcessingFunctional2D_NN<PRECOMP_T,PRECOMP_T>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field1,
        MultiNTensorField2D<PRECOMP_T>& field2,
        MultiNTensorField2D<int>& mask,
        plint level );

template
void applyProcessingFunctional<PRECOMP_T> (
        NTensorFieldBoxProcessingFunctional2D<PRECOMP_T>* functional,
        Box2D domain, std::vector<MultiNTensorField2D<PRECOMP_T>*> fields );
template
void integrateProcessingFunctional<PRECOMP_T> (
        NTensorFieldBoxProcessingFunctional2D<PRECOMP_T>* functional,
        Box2D domain, std::vector<MultiNTensorField2D<PRECOMP_T>*> fields, plint level );

template
void applyProcessingFunctional<PRECOMP_T> (
        ScalarFieldBoxProcessingFunctional2D<PRECOMP_T>* functional,
        Box2D domain, std::vector<MultiScalarField2D<PRECOMP_T>*> fields );
template
void integrateProcessingFunctional<PRECOMP_T> (
        ScalarFieldBoxProcessingFunctional2D<PRECOMP_T>* functional,
        Box2D domain, std::vector<MultiScalarField2D<PRECOMP_T>*> fields, plint level );

template
void applyProcessingFunctional<PRECOMP_T> (
        MaskedNTensorFieldBoxProcessingFunctional2D<PRECOMP_T>* functional,
        Box2D domain,
        std::vector<MultiNTensorField2D<PRECOMP_T>*> fields,
        MultiNTensorField2D<int>& mask);
template
void integrateProcessingFunctional<PRECOMP_T> (
        MaskedNTensorFieldBoxProcessingFunctional2D<PRECOMP_T>* functional,
        Box2D domain, std::vector<MultiNTensorField2D<PRECOMP_T>*> fields,
        MultiNTensorField2D<int>& mask, plint level );


/* *************** Bounded Boxed Data Processor functionals ****************** */

template void applyProcessingFunctional<PRECOMP_T> (
        BoundedBoxProcessingFunctional2D_N<PRECOMP_T>* functional,
        Box2D domain, MultiNTensorField2D<PRECOMP_T>& field, plint boundaryWidth );
template
void integrateProcessingFunctional<PRECOMP_T> (
        BoundedBoxProcessingFunctional2D_N<PRECOMP_T>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field, plint boundaryWidth, plint level );

template
void applyProcessingFunctional<PRECOMP_T> (
        BoundedMaskedBoxProcessingFunctional2D_N<PRECOMP_T>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field,
        MultiNTensorField2D<int>& mask,
        plint boundaryWidth );
template
void integrateProcessingFunctional<PRECOMP_T> (
        BoundedMaskedBoxProcessingFunctional2D_N<PRECOMP_T>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field,
        MultiNTensorField2D<int>& mask,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<PRECOMP_T,PRECOMP_T> (
        BoundedBoxProcessingFunctional2D_NN<PRECOMP_T,PRECOMP_T>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field1,
        MultiNTensorField2D<PRECOMP_T>& field2,
        plint boundaryWidth );
template
void integrateProcessingFunctional<PRECOMP_T,PRECOMP_T> (
        BoundedBoxProcessingFunctional2D_NN<PRECOMP_T,PRECOMP_T>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field1,
        MultiNTensorField2D<PRECOMP_T>& field2,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<PRECOMP_T,PRECOMP_T> (
        BoundedMaskedBoxProcessingFunctional2D_NN<PRECOMP_T,PRECOMP_T>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field1,
        MultiNTensorField2D<PRECOMP_T>& field2,
        MultiNTensorField2D<int>& mask,
        plint boundaryWidth );
template
void integrateProcessingFunctional<PRECOMP_T,PRECOMP_T> (
        BoundedMaskedBoxProcessingFunctional2D_NN<PRECOMP_T,PRECOMP_T>* functional,
        Box2D domain,
        MultiNTensorField2D<PRECOMP_T>& field1,
        MultiNTensorField2D<PRECOMP_T>& field2,
        MultiNTensorField2D<int>& mask,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<PRECOMP_T> (
        BoundedNTensorFieldBoxProcessingFunctional2D<PRECOMP_T>* functional,
        Box2D domain, std::vector<MultiNTensorField2D<PRECOMP_T>*> fields, plint boundaryWidth );
template
void integrateProcessingFunctional<PRECOMP_T> (
        BoundedNTensorFieldBoxProcessingFunctional2D<PRECOMP_T>* functional,
        Box2D domain, std::vector<MultiNTensorField2D<PRECOMP_T>*> fields, plint boundaryWidth, plint level );

}  // namespace plb
