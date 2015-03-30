/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
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

/* Main author: Orestis Malaspinas
 */

/** \file
 * Collision terms with external force -- generic code.
 */
#ifndef EXTERNAL_FORCE_DYNAMICS_HH
#define EXTERNAL_FORCE_DYNAMICS_HH

#include "basicDynamics/isoThermalDynamics.h"
#include "core/cell.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"
#include "latticeBoltzmann/d3q13Templates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/externalForceTemplates.h"
#include "core/latticeStatistics.h"
#include <algorithm>
#include "core/dynamicsIdentifiers.h"
#include <limits>

namespace plb {

/* *************** Class ExternalForceDynamics *********************************************** */

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
ExternalForceDynamics<T,Descriptor>::ExternalForceDynamics(T omega_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
void ExternalForceDynamics<T,Descriptor>::computeVelocity (
                                        Cell<T,Descriptor> const& cell, 
                                        Array<T,Descriptor<T>::d>& u ) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d> force, j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar, j);
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    
    T invRho = Descriptor<T>::invRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
        u[iD] = j[iD]*invRho + force[iD]/(T)2;
    }
    
}


/* *************** Class NaiveExternalForceBGKdynamics ********************************** */

template<typename T, template<typename U> class Descriptor>
int NaiveExternalForceBGKdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,NaiveExternalForceBGKdynamics<T,Descriptor> >("BGK_ExternalForce_Naive");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
NaiveExternalForceBGKdynamics<T,Descriptor>::NaiveExternalForceBGKdynamics(T omega_ )
    : ExternalForceDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
NaiveExternalForceBGKdynamics<T,Descriptor>* NaiveExternalForceBGKdynamics<T,Descriptor>::clone() const {
    return new NaiveExternalForceBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int NaiveExternalForceBGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void NaiveExternalForceBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d>  j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    externalForceTemplates<T,Descriptor>::addNaiveForce(cell);
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T NaiveExternalForceBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}



/* *************** Class NaiveExternalForcePrecondBGKdynamics ********************************** */

template<typename T, template<typename U> class Descriptor>
int NaiveExternalForcePrecondBGKdynamics<T,Descriptor>::id =
    meta::registerTwoParamDynamics<T,Descriptor,NaiveExternalForcePrecondBGKdynamics<T,Descriptor> >("Precond_BGK_ExternalForce_Naive");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
NaiveExternalForcePrecondBGKdynamics<T,Descriptor>::NaiveExternalForcePrecondBGKdynamics(T omega_, T invGamma_ )
    : ExternalForceDynamics<T,Descriptor>(omega_),
      invGamma(invGamma_)
{ }

template<typename T, template<typename U> class Descriptor>
NaiveExternalForcePrecondBGKdynamics<T,Descriptor>* NaiveExternalForcePrecondBGKdynamics<T,Descriptor>::clone() const {
    return new NaiveExternalForcePrecondBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int NaiveExternalForcePrecondBGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void NaiveExternalForcePrecondBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d>  j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = dynamicsTemplates<T,Descriptor>::precond_bgk_ma2_collision(cell, rhoBar, j, this->getOmega(), invGamma);
    externalForceTemplates<T,Descriptor>::addNaiveForce(cell);
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T NaiveExternalForcePrecondBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::precond_bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr, invGamma);
}


template<typename T, template<typename U> class Descriptor>
void NaiveExternalForcePrecondBGKdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(invGamma);
    serializer.addValue(this->getOmega());
}

template<typename T, template<typename U> class Descriptor>
void NaiveExternalForcePrecondBGKdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    invGamma = unserializer.readValue<T>();
    this->setOmega(unserializer.readValue<T>());
}


/* *************** Class GuoExternalForceBGKdynamics ********************************** */

template<typename T, template<typename U> class Descriptor>
int GuoExternalForceBGKdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,GuoExternalForceBGKdynamics<T,Descriptor> >("BGK_ExternalForce_Guo");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
GuoExternalForceBGKdynamics<T,Descriptor>::GuoExternalForceBGKdynamics(T omega_ )
    : ExternalForceDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceBGKdynamics<T,Descriptor>* GuoExternalForceBGKdynamics<T,Descriptor>::clone() const {
    return new GuoExternalForceBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int GuoExternalForceBGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar = this->computeRhoBar(cell);
    Array<T,Descriptor<T>::d> u, j;
    this->computeVelocity(cell, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
        j[iD] = rho * u[iD];
    }
    
    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    externalForceTemplates<T,Descriptor>::addGuoForce(cell, u, this->getOmega(), (T)1);
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T GuoExternalForceBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}


/* *************** Class ShanChenExternalForceBGKdynamics ********************************** */

template<typename T, template<typename U> class Descriptor>
int ShanChenExternalForceBGKdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,ShanChenExternalForceBGKdynamics<T,Descriptor> > 
                                      ("BGK_ExternalForce_ShanChen");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
ShanChenExternalForceBGKdynamics<T,Descriptor>::ShanChenExternalForceBGKdynamics(T omega_ )
    : ExternalForceDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
ShanChenExternalForceBGKdynamics<T,Descriptor>* ShanChenExternalForceBGKdynamics<T,Descriptor>::clone() const {
    return new ShanChenExternalForceBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int ShanChenExternalForceBGKdynamics<T,Descriptor>::getId() const {
    return id;
}


template<typename T, template<typename U> class Descriptor>
void ShanChenExternalForceBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    
    T invOmega = 1./this->getOmega();
    Array<T,Descriptor<T>::d> force;
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    Array<T,Descriptor<T>::d> jCorrected(j+invOmega*Descriptor<T>::fullRho(rhoBar) * force);
    
    dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rhoBar, jCorrected, this->getOmega());
    if (cell.takesStatistics()) {
        T uSqr = T();
        T invRho = Descriptor<T>::invRho(rhoBar);
        for (int iD=0; iD<Descriptor<T>::d; ++iD) {
            uSqr += util::sqr(j[iD]*invRho + 0.5*force[iD]);
        }
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T ShanChenExternalForceBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar ) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}


/* *************** Class HeExternalForceBGKdynamics ********************************** */

template<typename T, template<typename U> class Descriptor>
int HeExternalForceBGKdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,HeExternalForceBGKdynamics<T,Descriptor> >("BGK_ExternalForce_He");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
HeExternalForceBGKdynamics<T,Descriptor>::HeExternalForceBGKdynamics(T omega_ )
    : ExternalForceDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
HeExternalForceBGKdynamics<T,Descriptor>* HeExternalForceBGKdynamics<T,Descriptor>::clone() const {
    return new HeExternalForceBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int HeExternalForceBGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void HeExternalForceBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = externalForceTemplates<T,Descriptor>::heForcedBGKCollision (
            cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T HeExternalForceBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class IncGuoExternalForceBGKdynamics ********************************** */

template<typename T, template<typename U> class Descriptor>
int IncGuoExternalForceBGKdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,IncGuoExternalForceBGKdynamics<T,Descriptor> >("IncBGK_ExternalForce_Guo");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
IncGuoExternalForceBGKdynamics<T,Descriptor>::IncGuoExternalForceBGKdynamics(T omega_ )
    : ExternalForceDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
IncGuoExternalForceBGKdynamics<T,Descriptor>* IncGuoExternalForceBGKdynamics<T,Descriptor>::clone() const {
    return new IncGuoExternalForceBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int IncGuoExternalForceBGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void IncGuoExternalForceBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar = this->computeRhoBar(cell);
    Array<T,Descriptor<T>::d> u, j;
    this->computeVelocity(cell, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        j[iD] = rho*u[iD];
    }
    
    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_inc_collision(cell, rhoBar, j, this->getOmega());
    externalForceTemplates<T,Descriptor>::addGuoForce(cell, u, this->getOmega(), (T)1);
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T IncGuoExternalForceBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
bool IncGuoExternalForceBGKdynamics<T,Descriptor>::velIsJ() const {
    return true;
}

/* *************** Class ShanChenExternalForceRegularizedBGKdynamics ************************************ */

template<typename T, template<typename U> class Descriptor>
int ShanChenExternalForceRegularizedBGKdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,ShanChenExternalForceRegularizedBGKdynamics<T,Descriptor> >("Regularized_BGK_ExternalForce_Guo");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
ShanChenExternalForceRegularizedBGKdynamics<T,Descriptor>::ShanChenExternalForceRegularizedBGKdynamics(T omega_)
: ExternalForceDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
ShanChenExternalForceRegularizedBGKdynamics<T,Descriptor>* ShanChenExternalForceRegularizedBGKdynamics<T,Descriptor>::clone() const
{
    return new ShanChenExternalForceRegularizedBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int ShanChenExternalForceRegularizedBGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void ShanChenExternalForceRegularizedBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    
    T invOmega = 1./this->getOmega();
    Array<T,Descriptor<T>::d> force;
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    Array<T,Descriptor<T>::d> jCorrected(j+invOmega*Descriptor<T>::fullRho(rhoBar) * force);
    
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = dynamicsTemplates<T,Descriptor>::rlb_collision (
                 cell, rhoBar, invRho, jCorrected, PiNeq, this->getOmega() );
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T ShanChenExternalForceRegularizedBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}


}  // namespace plb

#endif  // EXTERNAL_FORCE_DYNAMICS_HH
