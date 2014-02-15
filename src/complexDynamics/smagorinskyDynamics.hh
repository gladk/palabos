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

/* Orestis Malaspinas designed some of the classes and concepts contained
 * in this file. */

#ifndef SMAGORINSKY_DYNAMICS_HH
#define SMAGORINSKY_DYNAMICS_HH

#include "complexDynamics/smagorinskyDynamics.h"
#include "core/util.h"
#include "core/latticeStatistics.h"
#include "core/dynamicsIdentifiers.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/externalForceTemplates.h"
#include "latticeBoltzmann/mrtTemplates.h"
#include <cmath>

namespace plb {

template<typename T, template<typename U> class Descriptor>
struct SmagoOperations {
    static T computePrefactor(T omega0, T cSmago) {
        return (T)0.5 * util::sqr(cSmago*omega0*Descriptor<T>::invCs2);
    }
    static T computeOmega(T omega0, T preFactor, T rhoBar, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq)
    {
        T PiNeqNormSqr = SymmetricTensor<T,Descriptor>::tensorNormSqr(PiNeq);
        T PiNeqNorm    = sqrt(PiNeqNormSqr);
        T alpha        = preFactor * Descriptor<T>::invRho(rhoBar);
        T linearTerm   = alpha*PiNeqNorm;
        T squareTerm   = (T)2*alpha*alpha*PiNeqNormSqr;
        // In the following formula, the square-root appearing in the explicit form of
        //   omega is developed to second-order.
        return omega0*(1-linearTerm+squareTerm);
    }
};

template<typename T, template<typename U> class Descriptor>
int SmagorinskyDynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,SmagorinskyDynamics<T,Descriptor> >("Smagorinsky_Generic");

template<typename T, template<typename U> class Descriptor>
SmagorinskyDynamics<T,Descriptor>::SmagorinskyDynamics(Dynamics<T,Descriptor>* baseDynamics_,
                                                       T omega0_, T cSmago_, bool automaticPrepareCollision)
    : OmegaFromPiDynamics<T,Descriptor>(baseDynamics_, automaticPrepareCollision),
      omega0(omega0_),
      cSmago(cSmago_),
      preFactor(SmagoOperations<T,Descriptor>::computePrefactor(omega0, cSmago))
{ }

template<typename T, template<typename U> class Descriptor>
SmagorinskyDynamics<T,Descriptor>::SmagorinskyDynamics(HierarchicUnserializer& unserializer)
    : OmegaFromPiDynamics<T,Descriptor>(0, false),
      omega0(T()),
      cSmago(T()),
      preFactor(T())
{
    unserialize(unserializer);
}


template<typename T, template<typename U> class Descriptor>
T SmagorinskyDynamics<T,Descriptor>::getOmegaFromPiAndRhoBar(Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T rhoBar) const
{
    return SmagoOperations<T,Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
}

template<typename T, template<typename U> class Descriptor>
SmagorinskyDynamics<T,Descriptor>* SmagorinskyDynamics<T,Descriptor>::clone() const {
    return new SmagorinskyDynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int SmagorinskyDynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyDynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
    OmegaFromPiDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyDynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
    preFactor = SmagoOperations<T,Descriptor>::computePrefactor(omega0, cSmago);
    OmegaFromPiDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyDynamics<T,Descriptor>::setOmega(T omega0_)
{
    // Just to be sure to avoid an undefined state, also reset the value of
    // omega in the baseDynamics, although the actual value of omega in the
    // baseDynamics is omega=omega0+deltaOmega and is recomputed before each
    // collision.
    omega0 = omega0_;
    preFactor = SmagoOperations<T,Descriptor>::computePrefactor(omega0, cSmago);
    this->getBaseDynamics().setOmega(omega0);
}

template<typename T, template<typename U> class Descriptor>
T SmagorinskyDynamics<T,Descriptor>::getOmega() const {
    return omega0;
}


/* *************** Class SmagorinskyBGKdynamics ************************************** */

template<typename T, template<typename U> class Descriptor>
int SmagorinskyBGKdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,SmagorinskyBGKdynamics<T,Descriptor> >("BGK_Smagorinsky");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
SmagorinskyBGKdynamics<T,Descriptor>::SmagorinskyBGKdynamics (
        T omega0_, T cSmago_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega0_),
      omega0(omega0_),
      cSmago(cSmago_),
      preFactor(SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago))
{ }

template<typename T, template<typename U> class Descriptor>
SmagorinskyBGKdynamics<T,Descriptor>::SmagorinskyBGKdynamics (
        HierarchicUnserializer& unserializer )
    : IsoThermalBulkDynamics<T,Descriptor>(T()),
      omega0(T()),
      cSmago(T()),
      preFactor(T())
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyBGKdynamics<T,Descriptor>::setOmega(T omega0_)
{
    omega0 = omega0_;
    preFactor = SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago);
}

template<typename T, template<typename U> class Descriptor>
T SmagorinskyBGKdynamics<T,Descriptor>::getOmega() const {
    return omega0;
}

template<typename T, template<typename U> class Descriptor>
SmagorinskyBGKdynamics<T,Descriptor>* SmagorinskyBGKdynamics<T,Descriptor>::clone() const {
    return new SmagorinskyBGKdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int SmagorinskyBGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyBGKdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyBGKdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
    preFactor = SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago);
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T,Descriptor>::computeOmega (
            omega0, preFactor, rhoBar, PiNeq );
    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rhoBar, j, omega);
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyBGKdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T thetaBar, BlockStatistics& stat)
{
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T,Descriptor>::computeOmega (
            omega0, preFactor, rhoBar, PiNeq );
    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rhoBar, j, omega);
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T SmagorinskyBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar ) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}


/* *************** Class GuoExternalForceSmagorinskyBGKdynamics ************************************** */

template<typename T, template<typename U> class Descriptor>
int GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor> >("BGK_Smagorinsky_Guo");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::GuoExternalForceSmagorinskyBGKdynamics (
        T omega0_, T cSmago_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega0_),
      omega0(omega0_),
      cSmago(cSmago_),
      preFactor(SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago))
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::GuoExternalForceSmagorinskyBGKdynamics (
        HierarchicUnserializer& unserializer )
    : IsoThermalBulkDynamics<T,Descriptor>(T()),
      omega0(T()),
      cSmago(T()),
      preFactor(T())
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::setOmega(T omega0_)
{
    omega0 = omega0_;
    preFactor = SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago);
}

template<typename T, template<typename U> class Descriptor>
T GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::getOmega() const {
    return omega0;
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::setParameter(plint whichParameter, T value)
{
    switch (whichParameter) {
        case dynamicParams::omega_shear:
            omega0 = value;
            break;
        case dynamicParams::smagorinskyConstant:
            cSmago = value;
            break;
        default:
            IsoThermalBulkDynamics<T,Descriptor>::setParameter(whichParameter, value);
    }

    preFactor = SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago);
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>* GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::clone() const {
    return new GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
    preFactor = SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago);
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::computeVelocity (
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

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T,Descriptor>::computeOmega (
            omega0, preFactor, rhoBar, PiNeq );
    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rhoBar, j, omega);

    Array<T,Descriptor<T>::d> force, u;
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    T invRho = Descriptor<T>::invRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        u[iD] = j[iD]*invRho + force[iD]/(T)2;
    }
    externalForceTemplates<T,Descriptor>::addGuoForce(cell, u, omega, (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T thetaBar, BlockStatistics& stat)
{
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T,Descriptor>::computeOmega (
            omega0, preFactor, rhoBar, PiNeq );
    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rhoBar, j, omega);

    Array<T,Descriptor<T>::d> force, u;
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    T invRho = Descriptor<T>::invRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        u[iD] = j[iD]*invRho + force[iD]/(T)2;
    }
    externalForceTemplates<T,Descriptor>::addGuoForce(cell, u, omega, (T)1);
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T GuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar ) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}



/* *************** Class SmagorinskyIncBGKdynamics ************************************** */


template<typename T, template<typename U> class Descriptor>
int SmagorinskyIncBGKdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,SmagorinskyIncBGKdynamics<T,Descriptor> >("IncBGK_Smagorinsky");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
SmagorinskyIncBGKdynamics<T,Descriptor>::SmagorinskyIncBGKdynamics (
        T omega0_, T cSmago_, T rho0 )
    : IsoThermalBulkDynamics<T,Descriptor>(omega0_),
      omega0(omega0_),
      cSmago(cSmago_),
      preFactor(SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago)),
      invRho0((T)1/rho0)
{ }

template<typename T, template<typename U> class Descriptor>
SmagorinskyIncBGKdynamics<T,Descriptor>::SmagorinskyIncBGKdynamics(HierarchicUnserializer& unserializer)
    : IsoThermalBulkDynamics<T,Descriptor>(T()),
      omega0(T()),
      cSmago(T()),
      preFactor(T()),
      invRho0(T())
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
SmagorinskyIncBGKdynamics<T,Descriptor>* SmagorinskyIncBGKdynamics<T,Descriptor>::clone() const {
    return new SmagorinskyIncBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int SmagorinskyIncBGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
    serializer.addValue(invRho0);
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
    unserializer.readValue(invRho0);
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
    preFactor = SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago);
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T,Descriptor>::setOmega(T omega0_)
{
    omega0 = omega0_;
}

template<typename T, template<typename U> class Descriptor>
T SmagorinskyIncBGKdynamics<T,Descriptor>::getOmega() const {
    return IsoThermalBulkDynamics<T,Descriptor>::getOmega();
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T,Descriptor>::setParameter(plint whichParameter, T value)
{
    switch(whichParameter) {
        case dynamicParams::omega_shear: omega0 = value; break;
        case 110: invRho0 = (T)1/value; break;
        case dynamicParams::smagorinskyConstant: cSmago = value; break;
        default: IsoThermalBulkDynamics<T,Descriptor>::setParameter(whichParameter, value);
    }
}

template<typename T, template<typename U> class Descriptor>
T SmagorinskyIncBGKdynamics<T,Descriptor>::getParameter(plint whichParameter) const {
    switch(whichParameter) {
        case dynamicParams::omega_shear: return omega0;
        case 110: return (T)1/invRho0;
        case dynamicParams::smagorinskyConstant: return cSmago;
        default: return IsoThermalBulkDynamics<T,Descriptor>::getParameter(whichParameter);
    }
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T,Descriptor>::computeVelocity (
                                        Cell<T,Descriptor> const& cell,
                                        Array<T,Descriptor<T>::d>& u ) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar, j);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
        u[iD] = invRho0*j[iD];
    }
}


template<typename T, template<typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T,Descriptor>::computeOmega (
            omega0, preFactor, rhoBar, PiNeq );
    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_inc_collision(cell, rhoBar, j, omega, invRho0);
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T thetaBar, BlockStatistics& stat)
{
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T,Descriptor>::computeOmega (
            omega0, preFactor, rhoBar, PiNeq );
    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_inc_collision(cell, rhoBar, j, omega, invRho0);
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T SmagorinskyIncBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar ) const
{
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho0, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
bool SmagorinskyIncBGKdynamics<T,Descriptor>::velIsJ() const {
    return true;
}


/* *************** Class IncGuoExternalForceSmagorinskyBGKdynamics ************************************** */


template<typename T, template<typename U> class Descriptor>
int IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor> >("IncBGK_Smagorinsky_Guo");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::IncGuoExternalForceSmagorinskyBGKdynamics (
        T omega0_, T cSmago_, T rho0 )
    : IsoThermalBulkDynamics<T,Descriptor>(omega0_),
      omega0(omega0_),
      cSmago(cSmago_),
      invRho0((T)1/rho0)
{ }

template<typename T, template<typename U> class Descriptor>
IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::IncGuoExternalForceSmagorinskyBGKdynamics(HierarchicUnserializer& unserializer)
    : IsoThermalBulkDynamics<T,Descriptor>(T()),
      omega0(T()),
      cSmago(T()),
      invRho0(T())
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>* IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::clone() const {
    return new IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
    serializer.addValue(invRho0);
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
    unserializer.readValue(invRho0);
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::setOmega(T omega0_)
{
    omega0 = omega0_;
}

template<typename T, template<typename U> class Descriptor>
T IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::getOmega() const {
    return IsoThermalBulkDynamics<T,Descriptor>::getOmega();
}

template<typename T, template<typename U> class Descriptor>
void IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::setParameter(plint whichParameter, T value)
{
    switch(whichParameter) {
        case dynamicParams::omega_shear: omega0 = value; break;
        case 110: invRho0 = (T)1/value; break;
        case dynamicParams::smagorinskyConstant: cSmago = value; break;
        default: IsoThermalBulkDynamics<T,Descriptor>::setParameter(whichParameter, value);
    }
}

template<typename T, template<typename U> class Descriptor>
T IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::getParameter(plint whichParameter) const {
    switch(whichParameter) {
        case dynamicParams::omega_shear: return omega0;
        case 110: return (T)1/invRho0;
        case dynamicParams::smagorinskyConstant: return cSmago;
        default: return IsoThermalBulkDynamics<T,Descriptor>::getParameter(whichParameter);
    }
}

template<typename T, template<typename U> class Descriptor>
void IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::computeVelocity (
                                        Cell<T,Descriptor> const& cell,
                                        Array<T,Descriptor<T>::d>& u ) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d> force, j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar, j);
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
        u[iD] = invRho0*(j[iD] + force[iD]/(T)2);
    }
}


template<typename T, template<typename U> class Descriptor>
void IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar = momentTemplates<T,Descriptor>::get_rhoBar(cell);
    Array<T,Descriptor<T>::d> u,j;
    computeVelocity(cell, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    j = u / invRho0;
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);

    T tau = (T)1/omega0;

    Array<T,Descriptor<T>::d> force; force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    Array<T,SymmetricTensor<T,Descriptor>::n> forceCorr; forceCorr.resetToZero();
    plint iPi = 0;
    for (plint iA = 0; iA < Descriptor<T>::d; ++iA) {
        for (plint iB = iA; iB < Descriptor<T>::d; ++iB) {
            forceCorr[iPi] = (u[iA]*force[iB]+u[iB]*force[iA]);
            ++iPi;
        }
    }
    PiNeq += forceCorr * (T)0.5;

    T normPiNeq = sqrt((T)2*SymmetricTensor<T,Descriptor>::tensorNormSqr(PiNeq));
    normPiNeq *= (T)2*rho*util::sqr(Descriptor<T>::cs2*cSmago*cSmago);

    T tauSGS = T();
    if (normPiNeq != T()) { // test to avoid division per 0
        Array<T,SymmetricTensor<T,Descriptor>::n> S =
                (-rho*tau*Descriptor<T>::cs2+sqrt(util::sqr(rho*tau*Descriptor<T>::cs2)+normPiNeq)) / normPiNeq * PiNeq;

        T sNorm = sqrt((T)2*SymmetricTensor<T,Descriptor>::tensorNormSqr(S));

        tauSGS = cSmago*cSmago*sNorm*Descriptor<T>::invCs2;
    }
    
    T omega = (T)1/(tau+tauSGS);

    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_inc_collision(cell, rhoBar, j, omega);
    
    externalForceTemplates<T,Descriptor>::addGuoForce(cell, u, omega, (T)1);
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T thetaBar, BlockStatistics& stat)
{
    // This function is not yet implemented. Actually, it is not even clear
    // what exactly the interpretation of this function should be: does j
    // include the f/2 correction of the body force or not?
    PLB_ASSERT( false );
}

template<typename T, template<typename U> class Descriptor>
T IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar ) const
{
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho0, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
bool IncGuoExternalForceSmagorinskyBGKdynamics<T,Descriptor>::velIsJ() const {
    return true;
}


/* *************** Class SmagorinskyRegularizedDynamics ************************************ */

template<typename T, template<typename U> class Descriptor>
int SmagorinskyRegularizedDynamics<T,Descriptor>::id =
    meta::registerTwoParamDynamics<T,Descriptor,SmagorinskyRegularizedDynamics<T,Descriptor> >("Smagorinsky_Regularized");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
SmagorinskyRegularizedDynamics<T,Descriptor>::SmagorinskyRegularizedDynamics (
        T omega0_, T cSmago_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega0_),
      omega0(omega0_),
      cSmago(cSmago_),
      preFactor(SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago))
{ }

template<typename T, template<typename U> class Descriptor>
void SmagorinskyRegularizedDynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
}


template<typename T, template<typename U> class Descriptor>
void SmagorinskyRegularizedDynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
    preFactor = SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago);
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyRegularizedDynamics<T,Descriptor>::setOmega(T omega0_)
{
    omega0 = omega0_;
    preFactor = SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago);
}

template<typename T, template<typename U> class Descriptor>
T SmagorinskyRegularizedDynamics<T,Descriptor>::getOmega() const {
    return omega0;
}

template<typename T, template<typename U> class Descriptor>
SmagorinskyRegularizedDynamics<T,Descriptor>* SmagorinskyRegularizedDynamics<T,Descriptor>::clone() const
{
    return new SmagorinskyRegularizedDynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int SmagorinskyRegularizedDynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyRegularizedDynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T,Descriptor>::computeOmega (
            omega0, preFactor, rhoBar, PiNeq );
    IsoThermalBulkDynamics<T,Descriptor>::setOmega(omega);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = dynamicsTemplates<T,Descriptor>::rlb_collision (
                 cell, rhoBar, invRho, j, PiNeq, omega );
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyRegularizedDynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar,
        Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat)
{
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T,Descriptor>::computeOmega (
            omega0, preFactor, rhoBar, PiNeq );
    IsoThermalBulkDynamics<T,Descriptor>::setOmega(omega);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = dynamicsTemplates<T,Descriptor>::rlb_collision (
                 cell, rhoBar, invRho, j, PiNeq, omega );
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T SmagorinskyRegularizedDynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}



/* *************** Class SecuredSmagorinskyRegularizedDynamics ************************************ */

template<typename T, template<typename U> class Descriptor>
int SecuredSmagorinskyRegularizedDynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,SecuredSmagorinskyRegularizedDynamics<T,Descriptor> >("Secured_Smagorinsky_Regularized");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
SecuredSmagorinskyRegularizedDynamics<T,Descriptor>::SecuredSmagorinskyRegularizedDynamics (
        T omega0_, T cSmago_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega0_),
      omega0(omega0_),
      cSmago(cSmago_),
      preFactor(SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago))
{ }
 
template<typename T, template<typename U> class Descriptor>
SecuredSmagorinskyRegularizedDynamics<T,Descriptor>::SecuredSmagorinskyRegularizedDynamics (
        HierarchicUnserializer& unserializer )
    : IsoThermalBulkDynamics<T,Descriptor>(T()),
      omega0(T()),
      cSmago(T()),
      preFactor(T())
{
    unserialize(unserializer);
}
    
template<typename T, template<typename U> class Descriptor>
void SecuredSmagorinskyRegularizedDynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
}


template<typename T, template<typename U> class Descriptor>
void SecuredSmagorinskyRegularizedDynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
    preFactor = SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago);
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void SecuredSmagorinskyRegularizedDynamics<T,Descriptor>::setOmega(T omega0_)
{
    omega0 = omega0_;
    preFactor = SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago);
}

template<typename T, template<typename U> class Descriptor>
T SecuredSmagorinskyRegularizedDynamics<T,Descriptor>::getOmega() const {
    return omega0;
}

template<typename T, template<typename U> class Descriptor>
SecuredSmagorinskyRegularizedDynamics<T,Descriptor>* SecuredSmagorinskyRegularizedDynamics<T,Descriptor>::clone() const
{
    return new SecuredSmagorinskyRegularizedDynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int SecuredSmagorinskyRegularizedDynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void SecuredSmagorinskyRegularizedDynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);

    T rho = Descriptor<T>::fullRho(rhoBar);
    T softLimit = 0.2;
    T hardLimit = 0.6;
    for (plint i=0; i<Descriptor<T>::d; ++i) {
        constrainValue(j[i], softLimit, hardLimit);
    }
    softLimit *= 0.1;
    hardLimit *= 0.1;
    for (plint i=0; i<SymmetricTensor<T,Descriptor>::n; ++i) {
        constrainValue(PiNeq[i], softLimit, hardLimit);
    }

    T omega = SmagoOperations<T,Descriptor>::computeOmega (
            omega0, preFactor, rhoBar, PiNeq );
    IsoThermalBulkDynamics<T,Descriptor>::setOmega(omega);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = dynamicsTemplates<T,Descriptor>::rlb_collision (
                 cell, rhoBar, invRho, j, PiNeq, omega );
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void SecuredSmagorinskyRegularizedDynamics<T,Descriptor>::constrainValue (
        T& value, T softLimit, T hardLimit )
{
    T fvalue = fabs(value);
    plint sign = fvalue > 0 ? +1 : -1;
    if (fvalue>softLimit) {
        if (fvalue>hardLimit) {
            value = softLimit*sign;
        }
        else {
            T diff = fvalue-softLimit;
            T correction = diff*diff / (hardLimit-softLimit);
            value -= correction*sign;
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void SecuredSmagorinskyRegularizedDynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar,
        Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat)
{
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);

    T rho = Descriptor<T>::fullRho(rhoBar);
    T softLimit = rho*0.1;
    T hardLimit = rho*0.2;
    Array<T,3> newJ(j);
    for (plint i=0; i<Descriptor<T>::d; ++i) {
        constrainValue(newJ[i], softLimit, hardLimit);
    }
    softLimit *= 0.05;
    hardLimit *= 0.05;
    for (plint i=0; i<SymmetricTensor<T,Descriptor>::n; ++i) {
        constrainValue(PiNeq[i], softLimit, hardLimit);
    }

    T omega = SmagoOperations<T,Descriptor>::computeOmega (
            omega0, preFactor, rhoBar, PiNeq );
    IsoThermalBulkDynamics<T,Descriptor>::setOmega(omega);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = dynamicsTemplates<T,Descriptor>::rlb_collision (
                 cell, rhoBar, invRho, newJ, PiNeq, omega );
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T SecuredSmagorinskyRegularizedDynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}


/* *************** Class SmagorinskyMRTdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
int SmagorinskyMRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,SmagorinskyMRTdynamics<T,Descriptor> >("MRT_Smagorinsky");

template<typename T, template<typename U> class Descriptor>
SmagorinskyMRTdynamics<T,Descriptor>::SmagorinskyMRTdynamics(plint externalParam_, T cSmago_)
    : IsoThermalBulkDynamics<T,Descriptor>(mrtParam<T,Descriptor>().get(externalParam_).getOmega()),
      param(0),
      externalParam(externalParam_),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
SmagorinskyMRTdynamics<T,Descriptor>::SmagorinskyMRTdynamics(MRTparam<T,Descriptor>* param_, T cSmago_)
    : IsoThermalBulkDynamics<T,Descriptor>(param_->getOmega()),
      param(param_),
      externalParam(-1),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
SmagorinskyMRTdynamics<T,Descriptor>::SmagorinskyMRTdynamics(HierarchicUnserializer& unserializer)
    : IsoThermalBulkDynamics<T,Descriptor>(T()),
      param(0),
      externalParam(-1)
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
SmagorinskyMRTdynamics<T,Descriptor>::SmagorinskyMRTdynamics(SmagorinskyMRTdynamics<T,Descriptor> const& rhs)
    : IsoThermalBulkDynamics<T,Descriptor>(rhs),
      param(rhs.param ? rhs.param->clone() : 0),
      externalParam(rhs.externalParam),
      cSmago(rhs.cSmago)
{ }

template<typename T, template<typename U> class Descriptor>
SmagorinskyMRTdynamics<T,Descriptor>&
    SmagorinskyMRTdynamics<T,Descriptor>::operator=(SmagorinskyMRTdynamics<T,Descriptor> const& rhs)
{
    SmagorinskyMRTdynamics<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
SmagorinskyMRTdynamics<T,Descriptor>::~SmagorinskyMRTdynamics()
{
    delete param;
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyMRTdynamics<T,Descriptor>::swap(SmagorinskyMRTdynamics<T,Descriptor>& rhs)
{
    std::swap(param, rhs.param);
    std::swap(externalParam, rhs.externalParam);
    T tmpOmega = this->getOmega();
    this->setOmega(rhs.getOmega());
    rhs.setOmega(tmpOmega);
    std::swap(cSmago,rhs.cSmago_);
}

template<typename T, template<typename U> class Descriptor>
SmagorinskyMRTdynamics<T,Descriptor>* SmagorinskyMRTdynamics<T,Descriptor>::clone() const {
    return new SmagorinskyMRTdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int SmagorinskyMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyMRTdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    int useExternalParamFlag = param ? 0:1;
    serializer.addValue(useExternalParamFlag);
    if (param) {
        param->serialize(serializer);
    }
    else {
        serializer.addValue(externalParam);
    }
    serializer.addValue(cSmago);
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyMRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    int useExternalParamFlag = unserializer.readValue<int>();
    if (useExternalParamFlag) {
        externalParam = unserializer.readValue<plint>();
        param = 0;
    }
    else {
        delete param;
        param = new MRTparam<T,Descriptor>(unserializer);
        externalParam = -1;
    }
    unserializer.readValue(cSmago);
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyMRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;

    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));

    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T rho = Descriptor<T>::fullRho(rhoBar);

    T tau = (T)1/parameter->getOmega();

    T normPiNeq = sqrt((T)2*SymmetricTensor<T,Descriptor>::tensorNormSqr(PiNeq));
    normPiNeq *= (T)2*rho*util::sqr(Descriptor<T>::cs2*cSmago*cSmago);

    T tauSGS = T();
    if (normPiNeq != T()) { // test to avoid division per 0
        Array<T,SymmetricTensor<T,Descriptor>::n> S =
                (-rho*tau*Descriptor<T>::cs2+sqrt(util::sqr(rho*tau*Descriptor<T>::cs2)+normPiNeq)) / normPiNeq * PiNeq;

        T sNorm = sqrt((T)2*SymmetricTensor<T,Descriptor>::tensorNormSqr(S));

        tauSGS = cSmago*cSmago*sNorm*Descriptor<T>::invCs2;
    }
    T omega = (T)1/(tau+tauSGS);

    std::vector<T> newS = parameter->getS();
    for (plint iPop  = 0; iPop < Descriptor<T>::shearIndexes; ++iPop) {
        newS[Descriptor<T>::shearViscIndexes[iPop]] = omega;
    }

    T invM_S[Descriptor<T>::q][Descriptor<T>::q];
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        for (plint jPop = 0; jPop < Descriptor<T>::q; ++jPop) {
            invM_S[iPop][jPop] = Descriptor<T>::invM[iPop][jPop] * newS[jPop];
        }
    }

    T jSqr = mrtTemp::mrtCollision(cell, rhoBar, j, invM_S);

    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr*invRho*invRho);
    }
}


template<typename T, template<typename U> class Descriptor>
void SmagorinskyMRTdynamics<T,Descriptor>::collideExternal (
            Cell<T,Descriptor>& cell, T rhoBar,
            Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat)
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;

    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);

    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));

    T rho = Descriptor<T>::fullRho(rhoBar);
    T tau = (T)1/parameter->getOmega();

    T normPiNeq = sqrt((T)2*SymmetricTensor<T,Descriptor>::tensorNormSqr(PiNeq));
    normPiNeq *= (T)2*rho*util::sqr(Descriptor<T>::cs2*cSmago*cSmago);

    T tauSGS = T();
    if (normPiNeq != T()) { // test to avoid division per 0
        Array<T,SymmetricTensor<T,Descriptor>::n> S =
                (-rho*tau*Descriptor<T>::cs2+sqrt(util::sqr(rho*tau*Descriptor<T>::cs2)+normPiNeq)) / normPiNeq * PiNeq;

        T sNorm = sqrt((T)2*SymmetricTensor<T,Descriptor>::tensorNormSqr(S));

        tauSGS = cSmago*cSmago*sNorm*Descriptor<T>::invCs2;
    }
    T omega = (T)1/(tau+tauSGS);

    std::vector<T> newS = parameter->getS();
    for (plint iPop  = 0; iPop < Descriptor<T>::shearIndexes; ++iPop) {
        newS[Descriptor<T>::shearViscIndexes[iPop]] = omega;
    }

    T invM_S[Descriptor<T>::q][Descriptor<T>::q];
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        for (plint jPop = 0; jPop < Descriptor<T>::q; ++jPop) {
            invM_S[iPop][jPop] = T();
            for (plint kPop = 0; kPop < Descriptor<T>::q; ++kPop) {
                if (kPop == jPop) {
                    invM_S[iPop][jPop] += Descriptor<T>::invM[iPop][kPop] * newS[kPop];
                }
            }
        }
    }

    T jSqr = mrtTemp::mrtCollision(cell, rhoBar, j, invM_S);

    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(stat, rhoBar, jSqr*invRho*invRho);
    }
}

template<typename T, template<typename U> class Descriptor>
T SmagorinskyMRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor> const& SmagorinskyMRTdynamics<T,Descriptor>::getMrtParameter() const {
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    return *parameter;
}

/* *************** Class ConsistentSmagorinskyMRTdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
int ConsistentSmagorinskyMRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,ConsistentSmagorinskyMRTdynamics<T,Descriptor> >("Consistent_MRT_Smagorinsky");

template<typename T, template<typename U> class Descriptor>
ConsistentSmagorinskyMRTdynamics<T,Descriptor>::ConsistentSmagorinskyMRTdynamics(plint externalParam_, T cSmago_)
    : IsoThermalBulkDynamics<T,Descriptor>(mrtParam<T,Descriptor>().get(externalParam_).getOmega()),
      param(0),
      externalParam(externalParam_),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
ConsistentSmagorinskyMRTdynamics<T,Descriptor>::ConsistentSmagorinskyMRTdynamics(MRTparam<T,Descriptor>* param_, T cSmago_)
    : IsoThermalBulkDynamics<T,Descriptor>(param_->getOmega()),
      param(param_),
      externalParam(-1),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
ConsistentSmagorinskyMRTdynamics<T,Descriptor>::ConsistentSmagorinskyMRTdynamics(HierarchicUnserializer& unserializer)
    : IsoThermalBulkDynamics<T,Descriptor>(T()),
      param(0),
      externalParam(-1)
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
ConsistentSmagorinskyMRTdynamics<T,Descriptor>::ConsistentSmagorinskyMRTdynamics(ConsistentSmagorinskyMRTdynamics<T,Descriptor> const& rhs)
    : IsoThermalBulkDynamics<T,Descriptor>(rhs),
      param(rhs.param ? rhs.param->clone() : 0),
      externalParam(rhs.externalParam),
      cSmago(rhs.cSmago)
{ }

template<typename T, template<typename U> class Descriptor>
ConsistentSmagorinskyMRTdynamics<T,Descriptor>&
    ConsistentSmagorinskyMRTdynamics<T,Descriptor>::operator=(ConsistentSmagorinskyMRTdynamics<T,Descriptor> const& rhs)
{
    ConsistentSmagorinskyMRTdynamics<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
ConsistentSmagorinskyMRTdynamics<T,Descriptor>::~ConsistentSmagorinskyMRTdynamics()
{
    delete param;
}

template<typename T, template<typename U> class Descriptor>
void ConsistentSmagorinskyMRTdynamics<T,Descriptor>::swap(ConsistentSmagorinskyMRTdynamics<T,Descriptor>& rhs)
{
    std::swap(param, rhs.param);
    std::swap(externalParam, rhs.externalParam);
    T tmpOmega = this->getOmega();
    this->setOmega(rhs.getOmega());
    rhs.setOmega(tmpOmega);
    std::swap(cSmago,rhs.cSmago_);
}

template<typename T, template<typename U> class Descriptor>
ConsistentSmagorinskyMRTdynamics<T,Descriptor>* ConsistentSmagorinskyMRTdynamics<T,Descriptor>::clone() const {
    return new ConsistentSmagorinskyMRTdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int ConsistentSmagorinskyMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void ConsistentSmagorinskyMRTdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    int useExternalParamFlag = param ? 0:1;
    serializer.addValue(useExternalParamFlag);
    if (param) {
        param->serialize(serializer);
    }
    else {
        serializer.addValue(externalParam);
    }
    serializer.addValue(cSmago);
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void ConsistentSmagorinskyMRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    int useExternalParamFlag = unserializer.readValue<int>();
    if (useExternalParamFlag) {
        externalParam = unserializer.readValue<plint>();
        param = 0;
    }
    else {
        delete param;
        param = new MRTparam<T,Descriptor>(unserializer);
        externalParam = -1;
    }
    unserializer.readValue(cSmago);
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void ConsistentSmagorinskyMRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    
    T rhoBar; 
    Array<T,Descriptor<T>::d> j;
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T rho = Descriptor<T>::fullRho(rhoBar);
    
    T tau = (T)1/parameter->getOmega();
    
    T normPiNeq = sqrt((T)2*SymmetricTensor<T,Descriptor>::tensorNormSqr(PiNeq));
    normPiNeq *= (T)2*rho*util::sqr(Descriptor<T>::cs2*cSmago*cSmago);
    
    Array<T,SymmetricTensor<T,Descriptor>::n> S; S.resetToZero();
    if (normPiNeq != T()) { // test to avoid division per 0
        S = (-rho*tau*Descriptor<T>::cs2+sqrt(util::sqr(rho*tau*Descriptor<T>::cs2)+normPiNeq)) / normPiNeq * PiNeq;
    }
    
    T jSqr = mrtTemp::smagorinskyMrtCollision(cell, rhoBar, j, S, cSmago, parameter->getInvM());
    
    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr*invRho*invRho);
    }
}


template<typename T, template<typename U> class Descriptor>
void ConsistentSmagorinskyMRTdynamics<T,Descriptor>::collideExternal (
            Cell<T,Descriptor>& cell, T rhoBar,
            Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat)
{
    PLB_ASSERT(false);
}

template<typename T, template<typename U> class Descriptor>
T ConsistentSmagorinskyMRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, (T)1, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor> const& ConsistentSmagorinskyMRTdynamics<T,Descriptor>::getMrtParameter() const {
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    return *parameter;
}


/* *************** Class ConsistentSmagorinskyQuasiIncMRTdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
int ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor> >("Consistent_QuasiIncMRT_Smagorinsky");

template<typename T, template<typename U> class Descriptor>
ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::ConsistentSmagorinskyQuasiIncMRTdynamics(plint externalParam_, T cSmago_)
    : IsoThermalBulkDynamics<T,Descriptor>(mrtParam<T,Descriptor>().get(externalParam_).getOmega()),
      param(0),
      externalParam(externalParam_),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::ConsistentSmagorinskyQuasiIncMRTdynamics(MRTparam<T,Descriptor>* param_, T cSmago_)
    : IsoThermalBulkDynamics<T,Descriptor>(param_->getOmega()),
      param(param_),
      externalParam(-1),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::ConsistentSmagorinskyQuasiIncMRTdynamics(HierarchicUnserializer& unserializer)
    : IsoThermalBulkDynamics<T,Descriptor>(T()),
      param(0),
      externalParam(-1)
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::ConsistentSmagorinskyQuasiIncMRTdynamics(ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor> const& rhs)
    : IsoThermalBulkDynamics<T,Descriptor>(rhs),
      param(rhs.param ? rhs.param->clone() : 0),
      externalParam(rhs.externalParam),
      cSmago(rhs.cSmago)
{ }

template<typename T, template<typename U> class Descriptor>
ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>&
    ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::operator=(ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor> const& rhs)
{
    ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::~ConsistentSmagorinskyQuasiIncMRTdynamics()
{
    delete param;
}

template<typename T, template<typename U> class Descriptor>
void ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::swap(ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>& rhs)
{
    std::swap(param, rhs.param);
    std::swap(externalParam, rhs.externalParam);
    T tmpOmega = this->getOmega();
    this->setOmega(rhs.getOmega());
    rhs.setOmega(tmpOmega);
    std::swap(cSmago,rhs.cSmago_);
}

template<typename T, template<typename U> class Descriptor>
ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>* ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::clone() const {
    return new ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    int useExternalParamFlag = param ? 0:1;
    serializer.addValue(useExternalParamFlag);
    if (param) {
        param->serialize(serializer);
    }
    else {
        serializer.addValue(externalParam);
    }
    serializer.addValue(cSmago);
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    int useExternalParamFlag = unserializer.readValue<int>();
    if (useExternalParamFlag) {
        externalParam = unserializer.readValue<plint>();
        param = 0;
    }
    else {
        delete param;
        param = new MRTparam<T,Descriptor>(unserializer);
        externalParam = -1;
    }
    unserializer.readValue(cSmago);
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    
    T rhoBar; 
    Array<T,Descriptor<T>::d> j;
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T rho = Descriptor<T>::fullRho(rhoBar);
    
    T tau = (T)1/parameter->getOmega();
    
    T normPiNeq = sqrt((T)2*SymmetricTensor<T,Descriptor>::tensorNormSqr(PiNeq));
    normPiNeq *= (T)2*rho*util::sqr(Descriptor<T>::cs2*cSmago*cSmago);
    
    Array<T,SymmetricTensor<T,Descriptor>::n> S; S.resetToZero();
    if (normPiNeq != T()) { // test to avoid division per 0
        S = (-rho*tau*Descriptor<T>::cs2+sqrt(util::sqr(rho*tau*Descriptor<T>::cs2)+normPiNeq)) / normPiNeq * PiNeq;
    }
    
    T jSqr = mrtTemp::quasiIncSmagorinskyMrtCollision(cell, rhoBar, j, S, cSmago, parameter->getInvM());
    
    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr*invRho*invRho);
    }
}


template<typename T, template<typename U> class Descriptor>
void ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::collideExternal (
            Cell<T,Descriptor>& cell, T rhoBar,
            Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat)
{
    PLB_ASSERT(false);
}

template<typename T, template<typename U> class Descriptor>
T ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, (T)1, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor> const& ConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::getMrtParameter() const {
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    return *parameter;
}

} // namespace plb

#endif  // SMAGORINSKY_DYNAMICS_HH

