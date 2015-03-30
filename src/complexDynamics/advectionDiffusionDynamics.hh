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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef ADVECTION_DIFFUSION_DYNAMICS_HH
#define ADVECTION_DIFFUSION_DYNAMICS_HH

#include "core/latticeStatistics.h"
#include "core/dynamicsIdentifiers.h"
#include "complexDynamics/advectionDiffusionDynamics.h"
#include "latticeBoltzmann/advectionDiffusionMomentTemplates.h"
#include "latticeBoltzmann/advectionDiffusionDynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/offEquilibriumAdvectionDiffusionTemplates.h"


namespace plb {

template<typename T, template<typename U> class Descriptor>
struct AD_SmagoOperations {
    static T computePrefactor(T omega0, T cSmago) {
        return util::sqr(cSmago*omega0*Descriptor<T>::invCs2);
    }
    static T computeOmega(T omega0, T alpha, Array<T,Descriptor<T>::d> const& j1)
    {
        T j1Norm      = norm(j1);
        T linearTerm  = alpha*j1Norm;
        T squareTerm  = linearTerm*linearTerm;
        // In the following formula, the square-root appearing in the explicit form of
        //   omega is developed to second-order.
        return omega0*(1-linearTerm+squareTerm);
    }
};

    
/* *************** Class AdvectionDiffusionDynamics ************************************ */

template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionDynamics<T,Descriptor>::AdvectionDiffusionDynamics(T omega_)
  : BasicBulkDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionDynamics<T,Descriptor>::regularize (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar ) const
{
    // jAdvDiff is the first order moment of 
    Array<T,Descriptor<T>::d> jEq;
    
    advectionDiffusionMomentTemplates<T,Descriptor>::get_jEq(cell, rhoBar, jEq);
    
    advectionDiffusionDynamicsTemplates<T,Descriptor>::regularize(cell,rhoBar,j,jEq);
}

template<typename T, template<typename U> class Descriptor>
T AdvectionDiffusionDynamics<T,Descriptor>::computeEbar(Cell<T,Descriptor> const& cell) const
{
    return T();
}


/* *************** Class SmagorinskyAdvectionDiffusionRLBdynamics *************** */

template<typename T, template<typename U> class Descriptor>
int SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor> >("SmagoAdvectionDiffusion_RLB");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>::SmagorinskyAdvectionDiffusionRLBdynamics (T omega_, T T0_, T cSmago_ )
    : AdvectionDiffusionDynamics<T,Descriptor>(omega_),
      invT0((T)1/T0_),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>::SmagorinskyAdvectionDiffusionRLBdynamics(HierarchicUnserializer& unserializer)
    : AdvectionDiffusionDynamics<T,Descriptor>((T)1)
{
    unserialize(unserializer);
}



template<typename T, template<typename U> class Descriptor>
SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>* SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>::clone() const {
    return new SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    // The order is important: it must be the same as the parameters of the
    // constructor, because otherwise the TwoParamGenerator fails.
    AdvectionDiffusionDynamics<T,Descriptor>::serialize(serializer);
    serializer.addValue(invT0);
    serializer.addValue(cSmago);
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    AdvectionDiffusionDynamics<T,Descriptor>::unserialize(unserializer);
    invT0 = unserializer.readValue<T>();
    cSmago = unserializer.readValue<T>();
}
 
template<typename T, template<typename U> class Descriptor>
void SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics ) 
{
    T rhoBar;
    Array<T,Descriptor<T>::d> jEq, jNeq;
    advectionDiffusionMomentTemplates<T,Descriptor>::get_rhoBar_jEq_jNeq(cell, rhoBar, jEq, jNeq);
    
    T omega0 = this->getOmega();
    T omega = omega0/(1.+(1.-omega0/2.0)*cSmago*cSmago*omega0*Descriptor<T>::invCs2*invT0*norm(jNeq));

    T uSqr = advectionDiffusionDynamicsTemplates<T,Descriptor>::
            no_corr_rlb_collision(cell, rhoBar, jEq, jNeq, omega );
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>::collideExternal (
            Cell<T,Descriptor>& cell, T rhoBar,
            Array<T,Descriptor<T>::d> const& jEq, T thetaBar, BlockStatistics& statistics )
{
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_j(cell, j);
    Array<T,Descriptor<T>::d> jNeq(j-jEq);

    T omega0 = this->getOmega();
    T omega = omega0/(1.+(1.-omega0/2.0)*cSmago*cSmago*omega0*Descriptor<T>::invCs2*invT0*norm(jNeq));

    T uSqr = advectionDiffusionDynamicsTemplates<T,Descriptor>::
            no_corr_rlb_collision(cell, rhoBar, jEq, jNeq, omega);
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}


/** \param j The parameter j is defined as j = j_advDiff = rho_advDiff*u_fluid
 */
template<typename T, template<typename U> class Descriptor>
T SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{ 
    return advectionDiffusionDynamicsTemplates<T,Descriptor>::
            bgk_ma1_equilibrium(iPop, rhoBar, j);
}


/* *************** Class AdvectionDiffusionRLBdynamics *************** */

template<typename T, template<typename U> class Descriptor>
int AdvectionDiffusionRLBdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,AdvectionDiffusionRLBdynamics<T,Descriptor> >("AdvectionDiffusion_RLB");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionRLBdynamics<T,Descriptor>::AdvectionDiffusionRLBdynamics (T omega_ )
    : AdvectionDiffusionDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionRLBdynamics<T,Descriptor>* AdvectionDiffusionRLBdynamics<T,Descriptor>::clone() const {
    return new AdvectionDiffusionRLBdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int AdvectionDiffusionRLBdynamics<T,Descriptor>::getId() const {
    return id;
}
 
template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionRLBdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics ) 
{
    T rhoBar;
    Array<T,Descriptor<T>::d> jEq, jNeq;
    advectionDiffusionMomentTemplates<T,Descriptor>::get_rhoBar_jEq_jNeq(cell, rhoBar, jEq, jNeq);
    
    T uSqr = advectionDiffusionDynamicsTemplates<T,Descriptor>::
            no_corr_rlb_collision(cell, rhoBar, jEq, jNeq, this->getOmega() );
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionRLBdynamics<T,Descriptor>::collideExternal (
            Cell<T,Descriptor>& cell, T rhoBar,
            Array<T,Descriptor<T>::d> const& jEq, T thetaBar, BlockStatistics& statistics )
{
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_j(cell, j);
    Array<T,Descriptor<T>::d> jNeq(j-jEq);

    T uSqr = advectionDiffusionDynamicsTemplates<T,Descriptor>::
            no_corr_rlb_collision(cell, rhoBar, jEq, jNeq, this->getOmega() );
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param j The parameter j is defined as j = j_advDiff = rho_advDiff*u_fluid
 */
template<typename T, template<typename U> class Descriptor>
T AdvectionDiffusionRLBdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{ 
    return advectionDiffusionDynamicsTemplates<T,Descriptor>::
            bgk_ma1_equilibrium(iPop, rhoBar, j);
}


/* *************** Class AdvectionDiffusionWithSourceRLBdynamics *************** */

template<typename T, template<typename U> class Descriptor>
int AdvectionDiffusionWithSourceRLBdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,AdvectionDiffusionWithSourceRLBdynamics<T,Descriptor> >("AdvectionDiffusionWithSource_RLB");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionWithSourceRLBdynamics<T,Descriptor>::AdvectionDiffusionWithSourceRLBdynamics (T omega_ )
    : AdvectionDiffusionDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionWithSourceRLBdynamics<T,Descriptor>* AdvectionDiffusionWithSourceRLBdynamics<T,Descriptor>::clone() const {
    return new AdvectionDiffusionWithSourceRLBdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int AdvectionDiffusionWithSourceRLBdynamics<T,Descriptor>::getId() const {
    return id;
}
 
template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionWithSourceRLBdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics ) 
{
    T rhoBar;
    Array<T,Descriptor<T>::d> jEq, jNeq;
    advectionDiffusionMomentTemplates<T,Descriptor>::get_rhoBar_jEq_jNeq(cell, rhoBar, jEq, jNeq);

    T sourceTerm = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);
    
    T uSqr = advectionDiffusionDynamicsTemplates<T,Descriptor>::
            no_corr_rlb_collision(cell, rhoBar, jEq, jNeq, this->getOmega(), sourceTerm );
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionWithSourceRLBdynamics<T,Descriptor>::collideExternal (
            Cell<T,Descriptor>& cell, T rhoBar,
            Array<T,Descriptor<T>::d> const& jEq, T thetaBar, BlockStatistics& statistics )
{
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_j(cell, j);
    Array<T,Descriptor<T>::d> jNeq(j-jEq);

    T sourceTerm = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);
    
    T uSqr = advectionDiffusionDynamicsTemplates<T,Descriptor>::
            no_corr_rlb_collision(cell, rhoBar, jEq, jNeq, this->getOmega(), sourceTerm );
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param j The parameter j is defined as j = j_advDiff = rho_advDiff*u_fluid
 */
template<typename T, template<typename U> class Descriptor>
T AdvectionDiffusionWithSourceRLBdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{ 
    return advectionDiffusionDynamicsTemplates<T,Descriptor>::
            bgk_ma1_equilibrium(iPop, rhoBar, j);
}


/* *************** Class AdvectionDiffusionBGKdynamics *************** */

template<typename T, template<typename U> class Descriptor>
int AdvectionDiffusionBGKdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,AdvectionDiffusionBGKdynamics<T,Descriptor> >("AdvectionDiffusion_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionBGKdynamics<T,Descriptor>::AdvectionDiffusionBGKdynamics (
        T omega_ )
    : AdvectionDiffusionDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionBGKdynamics<T,Descriptor>* AdvectionDiffusionBGKdynamics<T,Descriptor>::clone() const {
    return new AdvectionDiffusionBGKdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int AdvectionDiffusionBGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics ) 
{
    T rhoBar;
    Array<T,Descriptor<T>::d> jEq;
    advectionDiffusionMomentTemplates<T,Descriptor>::get_rhoBar_jEq(cell, rhoBar, jEq);
    
    T uSqr = advectionDiffusionDynamicsTemplates<T,Descriptor>::
            no_corr_bgk_collision(cell, rhoBar, jEq, this->getOmega());
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionBGKdynamics<T,Descriptor>::collideExternal (
            Cell<T,Descriptor>& cell, T rhoBar,
            Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& statistics )
{
    
    T uSqr = advectionDiffusionDynamicsTemplates<T,Descriptor>::
            no_corr_bgk_collision(cell, rhoBar, j, this->getOmega());
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param j The parameter j is defined as j = j_advDiff = rho_advDiff*u_fluid
 */
template<typename T, template<typename U> class Descriptor>
T AdvectionDiffusionBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{ 
    return advectionDiffusionDynamicsTemplates<T,Descriptor>::
            bgk_ma1_equilibrium(iPop, rhoBar, j);

}

// TODO implement decompose and recompose.
/* *************** Class CompleteAdvectionDiffusionBGKdynamics *************** */

template<typename T, template<typename U> class Descriptor>
int CompleteAdvectionDiffusionBGKdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,CompleteAdvectionDiffusionBGKdynamics<T,Descriptor> >("CompleteAdvectionDiffusion_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
CompleteAdvectionDiffusionBGKdynamics<T,Descriptor>::CompleteAdvectionDiffusionBGKdynamics (
        T omega_ )
    : AdvectionDiffusionDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
CompleteAdvectionDiffusionBGKdynamics<T,Descriptor>* CompleteAdvectionDiffusionBGKdynamics<T,Descriptor>::clone() const {
    return new CompleteAdvectionDiffusionBGKdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int CompleteAdvectionDiffusionBGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
T CompleteAdvectionDiffusionBGKdynamics<T,Descriptor>::computeDensity(Cell<T,Descriptor> const& cell) const
{
    T rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);
    return momentTemplates<T,Descriptor>::compute_rho(cell)*Descriptor<T>::invRho(rhoBar);
}

template<typename T, template<typename U> class Descriptor>
void CompleteAdvectionDiffusionBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics ) 
{
    T rhoPhiBar = momentTemplates<T,Descriptor>::get_rhoBar(cell);
    T rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);
    
    Array<T,Descriptor<T>::d> j;
    j.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
//     Array<T,SymmetricTensor<T,Descriptor>::n> piNeq;
//     piNeq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::piNeqBeginsAt));
    
    T rhoPhi = Descriptor<T>::fullRho(rhoPhiBar);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T phi    = rhoPhi * invRho;
    T invRhoPhiBar = Descriptor<T>::invRho(rhoPhiBar);
    
    T uSqr = dynamicsTemplates<T,Descriptor>::
            complete_bgk_ma2_collision(cell, rhoPhiBar, invRhoPhiBar, phi*j, this->getOmega());
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, Descriptor<T>::rhoBar(rhoPhi), uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void CompleteAdvectionDiffusionBGKdynamics<T,Descriptor>::collideExternal (
            Cell<T,Descriptor>& cell, T rhoBar,
            Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& statistics )
{
//     TODO IMPEMENT
    PLB_ASSERT(false);
}

/** \param j The parameter j is defined as j = j_advDiff = phi*rho_fluid*u_fluid
 */
template<typename T, template<typename U> class Descriptor>
T CompleteAdvectionDiffusionBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::complete_bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

// j = sum_i c_i*f_i, piNeq is nothing, jSqr also
template<typename T, template<typename U> class Descriptor>
void CompleteAdvectionDiffusionBGKdynamics<T,Descriptor>::regularize(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                            T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar) const
{
    T rhoPhiBar = rhoBar;
    rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);
    
    T phi = Descriptor<T>::fullRho(rhoPhiBar)*Descriptor<T>::invRho(rhoBar);

    Array<T,Descriptor<T>::d> jEq; 
    jEq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    
    jEq *= phi;
    
    Array<T,Descriptor<T>::d> jNeq = j - jEq;
    
    Array<T,SymmetricTensor<T,Descriptor>::n> pi;
    pi.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::piNeqBeginsAt)); 
    
    T omegaFluid = *cell.getExternal(Descriptor<T>::ExternalField::omegaBeginsAt);
    
    advectionDiffusionDynamicsTemplates<T,Descriptor>::complete_bgk_ma2_regularize(cell, rhoPhiBar, rhoBar, jEq, jNeq, pi, this->getOmega(), this->getOmega(), omegaFluid, omegaFluid );
} 


// TODO implement decompose and recompose.
/* *************** Class CompleteAdvectionDiffusionTRTdynamics *************** */

template<typename T, template<typename U> class Descriptor>
int CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::id =
    meta::registerTwoParamDynamics<T,Descriptor,CompleteAdvectionDiffusionTRTdynamics<T,Descriptor> >("CompleteAdvectionDiffusion_TRT");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::CompleteAdvectionDiffusionTRTdynamics (
        T omega_, T psi_ )
    : AdvectionDiffusionDynamics<T,Descriptor>(omega_), psi(psi_)
{ }

template<typename T, template<typename U> class Descriptor>
CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::CompleteAdvectionDiffusionTRTdynamics (
        T omega_)
    : AdvectionDiffusionDynamics<T,Descriptor>(omega_)
{ 
    psi = advectionDiffusionDynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor >::computePsiComplete(omega_);
}



template<typename T, template<typename U> class Descriptor>
CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>* CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::clone() const {
    return new CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    // The order is important: it must be the same as the parameters of the
    // constructor, because otherwise the TwoParamGenerator fails.
    AdvectionDiffusionDynamics<T,Descriptor>::serialize(serializer);
    serializer.addValue(psi);
}

template<typename T, template<typename U> class Descriptor>
void CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    AdvectionDiffusionDynamics<T,Descriptor>::unserialize(unserializer);
    psi = unserializer.readValue<T>();
}

template<typename T, template<typename U> class Descriptor>
T CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::computeDensity(Cell<T,Descriptor> const& cell) const
{
    T rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);
    return momentTemplates<T,Descriptor>::compute_rho(cell)*Descriptor<T>::invRho(rhoBar);
}

template<typename T, template<typename U> class Descriptor>
void CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics ) 
{
    T rhoPhiBar = momentTemplates<T,Descriptor>::get_rhoBar(cell);
    T rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);
    
    Array<T,Descriptor<T>::d> j;
    j.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    
    T rhoPhi = Descriptor<T>::fullRho(rhoPhiBar);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T phi    = rhoPhi * invRho;

    T uSqr = advectionDiffusionDynamicsTemplates<T,Descriptor>::complete_mrt_ma2_ext_rhoBar_j_collision(cell, rhoPhiBar, phi*j, this->getOmega(), psi);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, Descriptor<T>::rhoBar(rhoPhi), uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::collideExternal (
            Cell<T,Descriptor>& cell, T rhoBar,
            Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& statistics )
{
    T uSqr = advectionDiffusionDynamicsTemplates<T,Descriptor>::complete_mrt_ma2_ext_rhoBar_j_collision(cell, rhoBar, j, this->getOmega(), psi);
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param j The parameter j is defined as j = j_advDiff = phi*rho_fluid*u_fluid
 */
template<typename T, template<typename U> class Descriptor>
T CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::complete_bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

// j = sum_i c_i*f_i, piNeq is nothing, jSqr also
template<typename T, template<typename U> class Descriptor>
void CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::regularize(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                            T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar) const
{
    T rhoPhiBar = rhoBar;
    rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);
    
    T phi = Descriptor<T>::fullRho(rhoPhiBar)*Descriptor<T>::invRho(rhoBar);

    Array<T,Descriptor<T>::d> jEq; 
    jEq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    
    jEq *= phi;
    
    Array<T,Descriptor<T>::d> jNeq = j - jEq;
    
    Array<T,SymmetricTensor<T,Descriptor>::n> pi;
    pi.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::piNeqBeginsAt)); 
    
    T omegaFluid = *cell.getExternal(Descriptor<T>::ExternalField::omegaBeginsAt);
    
    advectionDiffusionDynamicsTemplates<T,Descriptor>::complete_bgk_ma2_regularize(cell, rhoPhiBar, rhoBar, jEq, jNeq, pi, this->getOmega(), psi, omegaFluid, omegaFluid );

} 


template<typename T, template<typename U> class Descriptor>
T CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::getParameter(plint whichParameter) const {
    switch (whichParameter) {
        case dynamicParams::omega_shear     : return this->getOmega();
        case dynamicParams::psi             : return getPsi();
    };
    return 0.;
}

template<typename T, template<typename U> class Descriptor>
void CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::setParameter(plint whichParameter, T value) {
    switch (whichParameter) {
        case dynamicParams::omega_shear     : this->setOmega(value);
        case dynamicParams::psi             : setPsi(value);
    };
}

template<typename T, template<typename U> class Descriptor>
T CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::getPsi() const {
    return psi;
}

template<typename T, template<typename U> class Descriptor>
void CompleteAdvectionDiffusionTRTdynamics<T,Descriptor>::setPsi(T psi_) {
    psi = psi_;
}

/* *************** Class AdvectionDiffusionWithSourceBGKdynamics *************** */

template<typename T, template<typename U> class Descriptor>
int AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor> >("AdvectionDiffusionWithSource_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::AdvectionDiffusionWithSourceBGKdynamics (
        T omega_ )
    : AdvectionDiffusionDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>* AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::clone() const {
    return new AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics ) 
{
    T rhoBar;
    Array<T,Descriptor<T>::d> jEq;
    advectionDiffusionMomentTemplates<T,Descriptor>::get_rhoBar_jEq(cell, rhoBar, jEq);

    T sourceTerm = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);
    
    T uSqr = advectionDiffusionDynamicsTemplates<T,Descriptor>::
            no_corr_bgk_collision(cell, rhoBar, jEq, this->getOmega(), sourceTerm);
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param j The parameter j is defined as j = j_advDiff = rho_advDiff*u_fluid
 */
template<typename T, template<typename U> class Descriptor>
T AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{ 
    return advectionDiffusionDynamicsTemplates<T,Descriptor>::
            bgk_ma1_equilibrium(iPop, rhoBar, j);

}


} // namespace plb

#endif  // ADVECTION_DIFFUSION_DYNAMICS_HH
