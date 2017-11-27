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
 * Interface to NumPy -- header file.
 */
#ifndef BLOCK_NUMPY_INTERFACE_2D_H
#define BLOCK_NUMPY_INTERFACE_2D_H

#include "multiBlock/multiDataField2D.h"

// http://blog.dhananjaynene.com/2009/03/constructor-method-overloading-in-python/

namespace plb {

template<typename T>
class NTensorField2NumPy2D {
public:
    NTensorField2NumPy2D(MultiNTensorField2D<T>& field_);
    NTensorField2NumPy2D(MultiNTensorField2D<T>& field_, Box2D const& domain_);
    void execute(T* array, int size);
    int getSize() const;
private:
    MultiNTensorField2D<T>& field;
    Box2D domain;
};

template<typename T>
class NumPy2NTensorField2D {
public:
    NumPy2NTensorField2D(MultiNTensorField2D<T>& field_);
    NumPy2NTensorField2D(MultiNTensorField2D<T>& field_, Box2D const& domain_);
    void execute(T* array, int size);
    int getSize() const;
private:
    MultiNTensorField2D<T>& field;
    Box2D domain;
};

}  // namespace plb

#endif  // NUMPY_INTERFACE_2D_H
