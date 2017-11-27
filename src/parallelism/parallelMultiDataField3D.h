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
 * Parallel access to elements of a scalar/tensor field -- header file.
 */
#ifndef PARALLEL_MULTI_DATA_FIELD_3D_H
#define PARALLEL_MULTI_DATA_FIELD_3D_H

#include "core/globalDefs.h"
#include "multiBlock/multiDataField3D.h"

#ifdef PLB_MPI_PARALLEL

namespace plb {

template<typename T>
class ParallelScalarAccess3D : public MultiScalarAccess3D<T> {
public:
    ParallelScalarAccess3D();
    virtual T& getDistributedScalar (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::map<plint,ScalarField3D<T>*>& fields );
    virtual T const& getDistributedScalar (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::map<plint,ScalarField3D<T>*> const& fields ) const;
    virtual ParallelScalarAccess3D<T>* clone() const;
private:
    mutable plint locatedBlock;
    mutable T distributedScalar;
};


template<typename T, int nDim>
class ParallelTensorAccess3D : public MultiTensorAccess3D<T,nDim> {
public:
    ParallelTensorAccess3D();
    virtual Array<T,nDim>& getDistributedTensor (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::map<plint,TensorField3D<T,nDim>*>& fields );
    virtual Array<T,nDim> const& getDistributedTensor (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::map<plint,TensorField3D<T,nDim>*> const& fields ) const;
    virtual ParallelTensorAccess3D<T,nDim>* clone() const;
private:
    mutable plint locatedBlock;
    mutable Array<T,nDim> distributedTensor;
};


template<typename T>
class ParallelNTensorAccess3D : public MultiNTensorAccess3D<T> {
public:
    ParallelNTensorAccess3D();
    virtual ~ParallelNTensorAccess3D();
    ParallelNTensorAccess3D(ParallelNTensorAccess3D<T> const& rhs);
    virtual T* getDistributedNTensor (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::map<plint,NTensorField3D<T>*>& fields );
    virtual T const* getDistributedNTensor (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::map<plint,NTensorField3D<T>*> const& fields ) const;
    virtual ParallelNTensorAccess3D<T>* clone() const;
private:
    ParallelNTensorAccess3D<T>& operator=(ParallelNTensorAccess3D<T> const& rhs) { return *this; }
private:
    mutable T* distributedNTensor;
    mutable plint locatedBlock;
};

}

#endif  // PLB_MPI_PARALLEL

#endif  // PARALLEL_MULTI_DATA_FIELD_3D_H
