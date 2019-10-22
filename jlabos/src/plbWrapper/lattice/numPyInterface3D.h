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
#ifndef NUMPY_INTERFACE_3D_H
#define NUMPY_INTERFACE_3D_H

#include "multiBlock/multiBlockLattice3D.h"

// http://blog.dhananjaynene.com/2009/03/constructor-method-overloading-in-python/

namespace plb {

template<typename T, template<typename U> class Descriptor>
class Lattice2NumPy3D {
public:
    Lattice2NumPy3D(MultiBlockLattice3D<T,Descriptor>& lattice_);
    Lattice2NumPy3D(MultiBlockLattice3D<T,Descriptor>& lattice_, Box3D const& domain_);
    void execute(T* array, int size);
    int getSize() const;
private:
    MultiBlockLattice3D<T,Descriptor>& lattice;
    Box3D domain;
};

template<typename T, template<typename U> class Descriptor>
class NumPy2Lattice3D {
public:
    NumPy2Lattice3D(MultiBlockLattice3D<T,Descriptor>& lattice_);
    NumPy2Lattice3D(MultiBlockLattice3D<T,Descriptor>& lattice_, Box3D const& domain_);
    void execute(T* array, int size);
    int getSize() const;
private:
    MultiBlockLattice3D<T,Descriptor>& lattice;
    Box3D domain;
};

}  // namespace plb

#endif  // NUMPY_INTERFACE_3D_H
