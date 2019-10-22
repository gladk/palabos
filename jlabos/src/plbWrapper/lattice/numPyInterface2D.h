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
#ifndef LATTICE_NUMPY_INTERFACE_2D_H
#define LATTICE_NUMPY_INTERFACE_2D_H

#include "multiBlock/multiBlockLattice2D.h"

// http://blog.dhananjaynene.com/2009/03/constructor-method-overloading-in-python/

namespace plb {

template<typename T, template<typename U> class Descriptor>
class Lattice2NumPy2D {
public:
    Lattice2NumPy2D(MultiBlockLattice2D<T,Descriptor>& lattice_);
    Lattice2NumPy2D(MultiBlockLattice2D<T,Descriptor>& lattice_, Box2D const& domain_);
    void execute(T* array, int size);
    int getSize() const;
private:
    MultiBlockLattice2D<T,Descriptor>& lattice;
    Box2D domain;
};

template<typename T, template<typename U> class Descriptor>
class NumPy2Lattice2D {
public:
    NumPy2Lattice2D(MultiBlockLattice2D<T,Descriptor>& lattice_);
    NumPy2Lattice2D(MultiBlockLattice2D<T,Descriptor>& lattice_, Box2D const& domain_);
    void execute(T* array, int size);
    int getSize() const;
private:
    MultiBlockLattice2D<T,Descriptor>& lattice;
    Box2D domain;
};

}  // namespace plb

#endif  // NUMPY_INTERFACE_2D_H
