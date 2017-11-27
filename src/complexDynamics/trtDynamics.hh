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

#ifndef TRT_DYNAMICS_HH
#define TRT_DYNAMICS_HH

#include "complexDynamics/trtDynamics.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "core/latticeStatistics.h"
#include <algorithm>
#include <limits>

namespace plb {

/* *************** Class TRTdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
const T TRTdynamics<T,Descriptor>::sMinus = 1.1;

template<typename T, template<typename U> class Descriptor>
int TRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,TRTdynamics<T,Descriptor> >("TRT");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
TRTdynamics<T,Descriptor>::TRTdynamics(T omega_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
TRTdynamics<T,Descriptor>::TRTdynamics(HierarchicUnserializer& unserializer)
    : IsoThermalBulkDynamics<T,Descriptor>(T())
{
    this->unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
TRTdynamics<T,Descriptor>* TRTdynamics<T,Descriptor>::clone() const {
    return new TRTdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int TRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void TRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    const T sPlus = this->getOmega();

    Array<T,Descriptor<T>::q> eq;
    // In the following, we number the plus/minus variables from 1 to (Q-1)/2.
    // So we allocate the index-zero memory location, and waste some memory
    // for convenience.
    Array<T,Descriptor<T>::q/2+1> eq_plus, eq_minus, f_plus, f_minus;

    Array<T,3> j;
    T rhoBar;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = normSqr(j);
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, eq);

    for (plint i=1; i<=Descriptor<T>::q/2; ++i) {
        eq_plus[i]  = 0.5*(eq[i] + eq[i+Descriptor<T>::q/2]);
        eq_minus[i] = 0.5*(eq[i] - eq[i+Descriptor<T>::q/2]);
        f_plus[i]   = 0.5*(cell[i] + cell[i+Descriptor<T>::q/2]);
        f_minus[i]  = 0.5*(cell[i] - cell[i+Descriptor<T>::q/2]);
    }

    cell[0] += -sPlus*cell[0] + sPlus*eq[0];

    for (plint i=1; i<=Descriptor<T>::q/2; ++i) {
        cell[i] += -sPlus*(f_plus[i]-eq_plus[i]) - sMinus*(f_minus[i]-eq_minus[i]);
        cell[i+Descriptor<T>::q/2] += -sPlus*(f_plus[i]-eq_plus[i]) + sMinus*(f_minus[i]-eq_minus[i]);
    }
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho );
    }
}

template<typename T, template<typename U> class Descriptor>
void TRTdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T thetaBar, BlockStatistics& statistics )
{
    const T sPlus = this->getOmega();

    Array<T,Descriptor<T>::q> eq;
    // In the following, we number the plus/minus variables from 1 to (Q-1)/2.
    // So we allocate the index-zero memory location, and waste some memory
    // for convenience.
    Array<T,Descriptor<T>::q/2+1> eq_plus, eq_minus, f_plus, f_minus;

    T jSqr = normSqr(j);
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, eq);

    for (plint i=1; i<=Descriptor<T>::q/2; ++i) {
        eq_plus[i]  = 0.5*(eq[i] + eq[i+Descriptor<T>::q/2]);
        eq_minus[i] = 0.5*(eq[i] - eq[i+Descriptor<T>::q/2]);
        f_plus[i]   = 0.5*(cell[i] + cell[i+Descriptor<T>::q/2]);
        f_minus[i]  = 0.5*(cell[i] - cell[i+Descriptor<T>::q/2]);
    }

    cell[0] += -sPlus*cell[0] + sPlus*eq[0];

    for (plint i=1; i<=Descriptor<T>::q/2; ++i) {
        cell[i] += -sPlus*(f_plus[i]-eq_plus[i]) - sMinus*(f_minus[i]-eq_minus[i]);
        cell[i+Descriptor<T>::q/2] += -sPlus*(f_plus[i]-eq_plus[i]) + sMinus*(f_minus[i]-eq_minus[i]);
    }
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar) );
    }
}

template<typename T, template<typename U> class Descriptor>
T TRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class IncTRTdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
const T IncTRTdynamics<T,Descriptor>::sMinus = 1.1;

template<typename T, template<typename U> class Descriptor>
int IncTRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,IncTRTdynamics<T,Descriptor> >("IncTRT");
    
/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
IncTRTdynamics<T,Descriptor>::IncTRTdynamics(T omega_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
IncTRTdynamics<T,Descriptor>::IncTRTdynamics(HierarchicUnserializer& unserializer)
    : IsoThermalBulkDynamics<T,Descriptor>(T())
{
    this->unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
IncTRTdynamics<T,Descriptor>* IncTRTdynamics<T,Descriptor>::clone() const {
    return new IncTRTdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int IncTRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void IncTRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    const T sPlus = this->getOmega();

    Array<T,Descriptor<T>::q> eq;
    // In the following, we number the plus/minus variables from 1 to (Q-1)/2.
    // So we allocate the index-zero memory location, and waste some memory
    // for convenience.
    Array<T,Descriptor<T>::q/2+1> eq_plus, eq_minus, f_plus, f_minus;

    Array<T,3> j;
    T rhoBar;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = normSqr(j);
    T invRho0 = 1.;
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria(rhoBar, invRho0, j, jSqr, eq);

    for (plint i=1; i<=Descriptor<T>::q/2; ++i) {
        eq_plus[i]  = 0.5*(eq[i] + eq[i+Descriptor<T>::q/2]);
        eq_minus[i] = 0.5*(eq[i] - eq[i+Descriptor<T>::q/2]);
        f_plus[i]   = 0.5*(cell[i] + cell[i+Descriptor<T>::q/2]);
        f_minus[i]  = 0.5*(cell[i] - cell[i+Descriptor<T>::q/2]);
    }

    cell[0] += -sPlus*cell[0] + sPlus*eq[0];

    for (plint i=1; i<=Descriptor<T>::q/2; ++i) {
        cell[i] += -sPlus*(f_plus[i]-eq_plus[i]) - sMinus*(f_minus[i]-eq_minus[i]);
        cell[i+Descriptor<T>::q/2] += -sPlus*(f_plus[i]-eq_plus[i]) + sMinus*(f_minus[i]-eq_minus[i]);
    }
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void IncTRTdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T thetaBar, BlockStatistics& statistics )
{
    const T sPlus = this->getOmega();

    Array<T,Descriptor<T>::q> eq;
    // In the following, we number the plus/minus variables from 1 to (Q-1)/2.
    // So we allocate the index-zero memory location, and waste some memory
    // for convenience.
    Array<T,Descriptor<T>::q/2+1> eq_plus, eq_minus, f_plus, f_minus;

    T jSqr = normSqr(j);
    T invRho0 = 1.;
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria(rhoBar, invRho0, j, jSqr, eq);

    for (plint i=1; i<=Descriptor<T>::q/2; ++i) {
        eq_plus[i]  = 0.5*(eq[i] + eq[i+Descriptor<T>::q/2]);
        eq_minus[i] = 0.5*(eq[i] - eq[i+Descriptor<T>::q/2]);
        f_plus[i]   = 0.5*(cell[i] + cell[i+Descriptor<T>::q/2]);
        f_minus[i]  = 0.5*(cell[i] - cell[i+Descriptor<T>::q/2]);
    }

    cell[0] += -sPlus*cell[0] + sPlus*eq[0];

    for (plint i=1; i<=Descriptor<T>::q/2; ++i) {
        cell[i] += -sPlus*(f_plus[i]-eq_plus[i]) - sMinus*(f_minus[i]-eq_minus[i]);
        cell[i+Descriptor<T>::q/2] += -sPlus*(f_plus[i]-eq_plus[i]) + sMinus*(f_minus[i]-eq_minus[i]);
    }
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T IncTRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, (T)1, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
bool IncTRTdynamics<T,Descriptor>::velIsJ() const {
    return true;
}

template<typename T, template<typename U> class Descriptor>
void IncTRTdynamics<T,Descriptor>::computeVelocity( Cell<T,Descriptor> const& cell,
                                                                 Array<T,Descriptor<T>::d>& u ) const 
{
    T dummyRhoBar;
    this->computeRhoBarJ(cell, dummyRhoBar, u);
}


}

#endif  // TRT_DYNAMICS_HH

