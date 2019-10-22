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
 * Coupling between grids of different refinement level -- header file.
 */

#ifndef MULTI_LEVEL_FIELD_GENERATOR_3D_H
#define MULTI_LEVEL_FIELD_GENERATOR_3D_H

#include "core/globalDefs.h"

#include "gridRefinement/octreeGridStructure.h"
#include "gridRefinement/couplingInterfaceGenerator3D.h"
#include "gridRefinement/multiLevelScalarField3D.h"
#include "gridRefinement/multiLevelTensorField3D.h"
#include "gridRefinement/multiLevelNTensorField3D.h"

namespace plb {

// ======================================================= //
// =================== Standard MultiLevels ============== //
// ======================================================= //

// ========== MultiLevelScalarField3D ==================== //
template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> > generateMultiLevelScalarField3D(
    const OctreeGridStructure &ogs);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, typename T3>
std::auto_ptr<MultiLevelScalarField3D<T3> > generateMultiLevelScalarField3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices);

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> > generateMultiLevelScalarField3D(
    const OctreeGridStructure &ogs, const Box3D &domain, plint levelOfDomain);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, 
    typename T3>
std::auto_ptr<MultiLevelScalarField3D<T3> > generateMultiLevelScalarField3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices, const Box3D &domain, plint levelOfDomain);

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> > generateMultiLevelScalarField3D(
    const OctreeGridStructure &ogs, T iniVal);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, typename T3>
std::auto_ptr<MultiLevelScalarField3D<T3> > generateMultiLevelScalarField3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices, T iniVal);

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> > generateMultiLevelScalarField3D(
    const OctreeGridStructure &ogs, const Box3D &domain, plint levelOfDomain, T iniVal);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, 
    typename T3>
std::auto_ptr<MultiLevelScalarField3D<T3> > generateMultiLevelScalarField3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices, const Box3D &domain, plint levelOfDomain, T iniVal);


// ========== MultiLevelTensorField3D ==================== //
template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> > generateMultiLevelTensorField3D(
    const OctreeGridStructure &ogs);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, 
    typename T3, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T3,nDim> > generateMultiLevelTensorField3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices);

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> > generateMultiLevelTensorField3D(
    const OctreeGridStructure &ogs, const Box3D &domain, plint levelOfDomain);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, 
    typename T3, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T3,nDim> > generateMultiLevelTensorField3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices, const Box3D &domain, plint levelOfDomain);

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> > generateMultiLevelTensorField3D(
    const OctreeGridStructure &ogs, const Array<T,nDim> &iniVal);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, 
    typename T3, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T3,nDim> > generateMultiLevelTensorField3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices, const Array<T,nDim> &iniVal);

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> > generateMultiLevelTensorField3D(
    const OctreeGridStructure &ogs, const Box3D &domain, plint levelOfDomain, const Array<T,nDim> &iniVal);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, 
    typename T3, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T3,nDim> > generateMultiLevelTensorField3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices, const Box3D &domain, plint levelOfDomain, const Array<T,nDim> &iniVal);

// ========== MultiLevelNTensorField3D ==================== //
template<typename T>
std::auto_ptr<MultiLevelNTensorField3D<T> > generateMultiLevelNTensorField3D(
    const OctreeGridStructure &ogs);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, 
    typename T3>
std::auto_ptr<MultiLevelNTensorField3D<T3> > generateMultiLevelNTensorField3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices);

template<typename T>
std::auto_ptr<MultiLevelNTensorField3D<T> > generateMultiLevelNTensorField3D(
    const OctreeGridStructure &ogs, const Box3D &domain, plint levelOfDomain);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, 
    typename T3>
std::auto_ptr<MultiLevelNTensorField3D<T3> > generateMultiLevelNTensorField3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices, const Box3D &domain, plint levelOfDomain);

// MultiLevels For Output
// ======================================================= //
// =================== Output MultiLevels ============== //
// ======================================================= //

// ========== MultiLevelScalarFieldForOutput3D ==================== //
template<typename T>
std::auto_ptr<MultiLevelScalarFieldForOutput3D<T> > generateMultiLevelScalarFieldForOutput3D(
    const OctreeGridStructure &ogs, bool crop);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, typename T3>
std::auto_ptr<MultiLevelScalarFieldForOutput3D<T3> > generateMultiLevelScalarFieldForOutput3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices, bool crop);

template<typename T>
std::auto_ptr<MultiLevelScalarFieldForOutput3D<T> > generateMultiLevelScalarFieldForOutput3D(
    const OctreeGridStructure &ogs, const Box3D &domain, plint levelOfDomain, bool crop);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, 
    typename T3>
std::auto_ptr<MultiLevelScalarFieldForOutput3D<T3> > generateMultiLevelScalarFieldForOutput3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices, const Box3D &domain, plint levelOfDomain, bool crop);


// ========== MultiLevelTensorFieldForOutput3D ==================== //
template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorFieldForOutput3D<T,nDim> > generateMultiLevelTensorFieldForOutput3D(
    const OctreeGridStructure &ogs, bool crop);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, 
    typename T3, int nDim>
std::auto_ptr<MultiLevelTensorFieldForOutput3D<T3,nDim> > generateMultiLevelTensorFieldForOutput3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices, bool crop);

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorFieldForOutput3D<T,nDim> > generateMultiLevelTensorFieldForOutput3D(
    const OctreeGridStructure &ogs, const Box3D &domain, plint levelOfDomain, bool crop);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, 
    typename T3, int nDim>
std::auto_ptr<MultiLevelTensorFieldForOutput3D<T3,nDim> > generateMultiLevelTensorFieldForOutput3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices, const Box3D &domain, plint levelOfDomain, bool crop);

// ========== MultiLevelNTensorFieldForOutput3D ==================== //
template<typename T>
std::auto_ptr<MultiLevelNTensorFieldForOutput3D<T> > generateMultiLevelNTensorFieldForOutput3D(
    const OctreeGridStructure &ogs, bool crop);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, 
    typename T3>
std::auto_ptr<MultiLevelNTensorFieldForOutput3D<T3> > generateMultiLevelNTensorFieldForOutput3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices, bool crop);

template<typename T>
std::auto_ptr<MultiLevelNTensorFieldForOutput3D<T> > generateMultiLevelNTensorFieldForOutput3D(
    const OctreeGridStructure &ogs, const Box3D &domain, plint levelOfDomain, bool crop);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine, 
    typename T3>
std::auto_ptr<MultiLevelNTensorFieldForOutput3D<T3> > generateMultiLevelNTensorFieldForOutput3D(
    const MultiLevelCoupling3D<T,Descriptor,Engine> &lattices, const Box3D &domain, plint levelOfDomain, bool crop);

}  // namespace plb

#endif  // MULTI_LEVEL_FIELD_GENERATOR_3D_H
