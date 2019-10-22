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
 * Base class for scalar, vector and tensor fields for 2D data analysis -- header file.
 */

#ifndef DATA_FIELD_BASE_2D_H
#define DATA_FIELD_BASE_2D_H

#include "core/globalDefs.h"
#include "core/block2D.h"
#include "core/geometry2D.h"
#include "core/array.h"

namespace plb {

/// Interface for the variants of 2D scalar fields.
template<typename T>
class ScalarFieldBase2D {
public:
    virtual ~ScalarFieldBase2D() { }
public:
    virtual void reset() =0;
    virtual T& get(plint iX, plint iY) =0;
    virtual T const& get(plint iX, plint iY) const =0;
};

/// Interface for the variants of 2D vector and tensor fields.
template<typename T, int nDim>
class TensorFieldBase2D {
public:
    virtual ~TensorFieldBase2D() { }
public:
    virtual void reset() =0;
    virtual Array<T,nDim>& get(plint iX, plint iY) =0;
    virtual Array<T,nDim> const& get(plint iX, plint iY) const =0;
};

/// Interface for the variants of generic-sized 2D vector and tensor fields.
/** The main purpose for these classes is use in dynamically typed languages
 *  like Python. In C++ it's most often better to use the static-sized
 *  TensorField to guarantee type safety.
 */
template<typename T>
class NTensorFieldBase2D {
public:
    NTensorFieldBase2D(int ndim_) : ndim(ndim_) { }
    NTensorFieldBase2D(NTensorFieldBase2D<T> const& rhs) : ndim(rhs.ndim) { }
    void swap(NTensorFieldBase2D& rhs) { std::swap(ndim, rhs.ndim); }
    virtual ~NTensorFieldBase2D() { }
public:
    virtual void reset() =0;
    virtual T* get(plint iX, plint iY) =0;
    virtual T const* get(plint iX, plint iY) const =0;
    plint getNdim() const { return ndim; }
private:
    plint ndim;
};

}  // namespace plb

#endif  // DATA_FIELD_BASE_2D_H
