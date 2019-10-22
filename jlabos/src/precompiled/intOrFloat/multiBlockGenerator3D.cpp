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
 * template instantiation.
 */
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiBlock/multiBlockGenerator3D.hh"

namespace plb {

template
std::auto_ptr<MultiNTensorField3D<PRECOMP_T> > defaultGenerateMultiNTensorField3D<PRECOMP_T> (
        MultiBlockManagement3D const& management, plint nDim );

template
MultiNTensorField3D<PRECOMP_T>* generateMultiNTensorField3D(Box3D const& domain, plint ndim);

template
MultiNTensorField3D<int>*
    generateNTensorFieldFromNTensor3D (
            MultiNTensorField3D<PRECOMP_T> const& field,
            Box3D const& intersection, plint nDim);

template
MultiNTensorField3D<PRECOMP_T>* clone<PRECOMP_T> (
        MultiNTensorField3D<PRECOMP_T>& originalField, Box3D const& subDomain, bool crop );

template
MultiNTensorField3D<PRECOMP_T>* generateMultiNTensorField<PRECOMP_T> (
        MultiBlock3D const& originalField, Box3D const& intersection,
        plint nDim, bool crop );

template
MultiNTensorField3D<PRECOMP_T>* generateIntersectMultiNTensorField<PRECOMP_T> (
        MultiBlock3D const& originalField1,
        MultiBlock3D const& originalField2,
        plint nDim, bool crop );

template
MultiNTensorField3D<PRECOMP_T>* generateIntersectMultiNTensorField<PRECOMP_T> (
        MultiBlock3D const& originalField1,
        MultiBlock3D const& originalField2,
        Box3D const& intersection, plint nDim,  bool crop );

template
MultiNTensorField3D<PRECOMP_T>* generateJoinMultiNTensorField<PRECOMP_T> (
        MultiBlock3D const& originalField1,
        MultiBlock3D const& originalField2,
        plint nDim );

template
MultiNTensorField3D<PRECOMP_T>* extend<PRECOMP_T> (
        MultiNTensorField3D<PRECOMP_T>& originalBlock, Box3D const& addedBlock );

template
MultiNTensorField3D<PRECOMP_T>* except<PRECOMP_T> (
        MultiNTensorField3D<PRECOMP_T>& originalBlock,
        Box3D const& exceptedBlock );

template
MultiNTensorField3D<PRECOMP_T>* align<PRECOMP_T> (
        MultiNTensorField3D<PRECOMP_T> const& originalBlock,
        MultiBlock3D const& partnerBlock );

template
MultiNTensorField3D<PRECOMP_T>* reparallelize<PRECOMP_T> (
        MultiNTensorField3D<PRECOMP_T> const& originalBlock );


template
std::auto_ptr<MultiScalarField3D<PRECOMP_T> > defaultGenerateMultiScalarField3D<PRECOMP_T> (
        MultiBlockManagement3D const& management, PRECOMP_T iniVal );

}  // namespace plb
