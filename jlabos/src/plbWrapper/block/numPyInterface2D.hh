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
 * Interface to NumPy -- generic code.
 */
#ifndef NUMPY_INTERFACE_2D_HH
#define NUMPY_INTERFACE_2D_HH

#include "plbWrapper/block/numPyInterface2D.h"
#include "core/serializer.h"
#include "parallelism/mpiManager.h"
#include "io/parallelIO.h"

namespace plb {

/* *************** Class NTensorField2NumPy2D ************************************ */

template<typename T>
NTensorField2NumPy2D<T>::NTensorField2NumPy2D(MultiNTensorField2D<T>& field_)
    : field(field_),
      domain(field.getBoundingBox())
{ }

template<typename T>
NTensorField2NumPy2D<T>::NTensorField2NumPy2D(MultiNTensorField2D<T>& field_, Box2D const& domain_)
    : field(field_),
      domain(domain_)
{ }

template<typename T>
void NTensorField2NumPy2D<T>::execute(T* array, int size) {
    serializerToSink (
            field.getBlockSerializer(domain, IndexOrdering::forward),
            new WriteToSerialArray<T>(array, (pluint)size)
    );
    global::mpi().bCast(array, size);
}

template<typename T>
int NTensorField2NumPy2D<T>::getSize() const {
    return (int) (domain.nCells() * field.getNdim());
}


/* *************** Class NumPy2NTensorField2D ************************************ */

template<typename T>
NumPy2NTensorField2D<T>::NumPy2NTensorField2D(MultiNTensorField2D<T>& field_)
    : field(field_),
      domain(field.getBoundingBox())
{ }

template<typename T>
NumPy2NTensorField2D<T>::NumPy2NTensorField2D(MultiNTensorField2D<T>& field_, Box2D const& domain_)
    : field(field_),
      domain(domain_)
{ }

template<typename T>
void NumPy2NTensorField2D<T>::execute(T* array, int size) {
    sourceToUnSerializer (
            new ReadFromSerialArray<T>(array, (pluint)size),
            field.getBlockUnSerializer(domain, IndexOrdering::forward)
    );
}

template<typename T>
int NumPy2NTensorField2D<T>::getSize() const {
    return (int) (domain.nCells() * field.getNdim());
}

}  // namespace plb

#endif  // NUMPY_INTERFACE_2D_HH
