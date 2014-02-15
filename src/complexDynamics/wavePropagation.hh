/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2013 FlowKit Sarl
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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef WAVE_PROPAGATION_HH
#define WAVE_PROPAGATION_HH

#include "basicDynamics/isoThermalDynamics.h"
#include "core/cell.h"
#include "core/dynamicsIdentifiers.h"
#include "core/latticeStatistics.h"

namespace plb {

/* *************** Class WaveDynamics ************************************ */

template<typename T, template<typename U> class Descriptor>
int WaveDynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,WaveDynamics<T,Descriptor> >("Wave");

/** \param vs2_ speed of sound
 *  \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
WaveDynamics<T,Descriptor>::WaveDynamics(T vs2_)
    : IsoThermalBulkDynamics<T,Descriptor>((T)2.0),
      vs2(vs2_)
{ }

template<typename T, template<typename U> class Descriptor>
void WaveDynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(vs2);
}

template<typename T, template<typename U> class Descriptor>
void WaveDynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    vs2 = unserializer.readValue<T>();
}

template<typename T, template<typename U> class Descriptor>
WaveDynamics<T,Descriptor>* WaveDynamics<T,Descriptor>::clone() const
{
    return new WaveDynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int WaveDynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void WaveDynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = waveCollision(cell, rhoBar, j, vs2);
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void WaveDynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, T rhoBar,
        Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat )
{
    waveCollision(cell, rhoBar, j, vs2);
}

template<typename T, template<typename U> class Descriptor>
T WaveDynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return waveEquilibrium(iPop, rhoBar, invRho, j, jSqr, vs2);
}

template<typename T, template<typename U> class Descriptor>
T WaveDynamics<T,Descriptor>::getParameter(plint whichParameter) const {
    switch (whichParameter) {
        case dynamicParams::sqrSpeedOfSound : return this->getVs2();
    };
    return 0.;
}

template<typename T, template<typename U> class Descriptor>
void WaveDynamics<T,Descriptor>::setParameter(plint whichParameter, T value) {
    switch (whichParameter) {
        case dynamicParams::sqrSpeedOfSound : setVs2(value);
    };
}

template<typename T, template<typename U> class Descriptor>
T WaveDynamics<T,Descriptor>::getVs2() const {
    return vs2;
}

template<typename T, template<typename U> class Descriptor>
void WaveDynamics<T,Descriptor>::setVs2(T vs2_) {
    vs2 = vs2_;
}

template<typename T, template<typename U> class Descriptor>
T WaveDynamics<T,Descriptor>::waveCollision (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T vs2)
{
    const T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    T invRho = Descriptor<T>::invRho(rhoBar);
    for (plint iPop=0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] *= (T)-1.0;
        cell[iPop] += (T)2.0 * waveEquilibrium(iPop, rhoBar, invRho, j, jSqr, vs2);
    }
    return invRho*invRho*jSqr;
}

template<typename T, template<typename U> class Descriptor>
T WaveDynamics<T,Descriptor>::waveEquilibrium (
        plint iPop, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j, T jSqr, T vs2)
{
    T kappa = vs2 - Descriptor<T>::cs2;
    if (iPop==0) {
        return Descriptor<T>::invCs2 * (
                     kappa * (Descriptor<T>::t[0]-(T)1)
                   + rhoBar * (Descriptor<T>::t[0]*vs2-kappa)) ;
    }
    else {
        T c_j = T();
        for (int iD=0; iD < Descriptor<T>::d; ++iD) {
           c_j += Descriptor<T>::c[iPop][iD]*j[iD];
        }
        return Descriptor<T>::invCs2 * Descriptor<T>::t[iPop] * (
                     kappa + rhoBar * vs2 + c_j );
    }
}

}  // namespace plb

#endif  // WAVE_PROPAGATION_HH

