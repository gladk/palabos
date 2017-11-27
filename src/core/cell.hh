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
 * Definition of a LB cell -- generic implementation.
 */
#ifndef CELL_HH
#define CELL_HH

#include <algorithm>
#include "core/cell.h"
#include "core/util.h"
#include <cstring>

namespace plb {

////////////////////////// Class Cell /////////////////////////////

/** The possibility to default construct Cell objects facilitates
 * their use in various types of containers. However, they can not
 * be used directly after construction; the method attributeDynamics()
 * must be used first.
 */
template<typename T, template<typename U> class Descriptor>
Cell<T,Descriptor>::Cell()
    : takesStat(true), dynamics(0)
{
    iniPop();
    iniExternal();
}

/** This constructor initializes the dynamics, but not the values
 * of the distribution functions. Remember that the dynamics is not
 * owned by the Cell object, the user must ensure its proper
 * destruction and a sufficient life time.
 */
template<typename T, template<typename U> class Descriptor>
Cell<T,Descriptor>::Cell(Dynamics<T,Descriptor>* dynamics_)
    : takesStat(true), dynamics(dynamics_)
{
    iniPop();
    iniExternal();
}

template<typename T, template<typename U> class Descriptor>
void Cell<T,Descriptor>::attributeDynamics(Dynamics<T,Descriptor>* dynamics_) {
    dynamics = dynamics_;
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor> const& Cell<T,Descriptor>::getDynamics() const {
    PLB_PRECONDITION(dynamics);
    return *dynamics;
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>& Cell<T,Descriptor>::getDynamics() {
    PLB_PRECONDITION(dynamics);
    return *dynamics;
}

template<typename T, template<typename U> class Descriptor>
void Cell<T,Descriptor>::revert() {
    for (plint iPop=1; iPop<=Descriptor<T>::numPop/2; ++iPop) {
        std::swap(f[iPop],f[iPop+Descriptor<T>::numPop/2]);
    }
}

template<typename T, template<typename U> class Descriptor>
void Cell<T,Descriptor>::iniPop() {
    f.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void Cell<T,Descriptor>::iniExternal() {
    for (plint iData=0; iData<Descriptor<T>::ExternalField::numScalars; ++iData) {
        *external.get(iData) = T();
    }
}

template<typename T, template<typename U> class Descriptor>
void Cell<T,Descriptor>::serialize(char* data) const {
    const plint numPop = Descriptor<T>::numPop;
    const plint numExt = Descriptor<T>::ExternalField::numScalars;
    memcpy((void*)data, (const void*)(&f[0]), numPop*sizeof(T));
    if (numExt>0) {
        memcpy((void*)(data+numPop*sizeof(T)), (const void*)(external.get(0)), numExt*sizeof(T));
    }
}

template<typename T, template<typename U> class Descriptor>
void Cell<T,Descriptor>::unSerialize(char const* data) {
    const plint numPop = Descriptor<T>::numPop;
    const plint numExt = Descriptor<T>::ExternalField::numScalars;

    memcpy((void*)(&f[0]), (const void*)data, numPop*sizeof(T));
    if (numExt>0) {
        memcpy((void*)(external.get(0)), (const void*)(data+numPop*sizeof(T)), numExt*sizeof(T));
    }
}

template<typename T, template<typename U> class Descriptor>
void iniCellAtEquilibrium(Cell<T,Descriptor>& cell, T density, Array<T,Descriptor<T>::d> const& velocity) {
    Array<T,Descriptor<T>::d> j;
    VectorTemplate<T,Descriptor>::multiplyByScalar(velocity, density, j);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    T rhoBar = Descriptor<T>::rhoBar(density);
    for (plint iPop=0; iPop<Descriptor<T>::numPop; ++iPop) {
        cell[iPop] = cell.computeEquilibrium(iPop, rhoBar, j, jSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void iniCellAtEquilibrium(Cell<T,Descriptor>& cell, T density, Array<T,Descriptor<T>::d> const& velocity, T temperature) {
    Array<T,Descriptor<T>::d> j;
    VectorTemplate<T,Descriptor>::multiplyByScalar(velocity, density, j);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    T rhoBar = Descriptor<T>::rhoBar(density);
    for (plint iPop=0; iPop<Descriptor<T>::numPop; ++iPop) {
        cell[iPop] = cell.computeEquilibrium(iPop, rhoBar, j, jSqr, temperature-(T)1);
    }
}

}  // namespace plb

#endif  // CELL_HH
