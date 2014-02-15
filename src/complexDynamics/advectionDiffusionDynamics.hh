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
#include "latticeBoltzmann/offEquilibriumAdvectionDiffusionTemplates.h"


namespace plb {

template<typename T, template<typename U> class Descriptor>
struct AD_SmagoOperations {
    static T computePrefactor(T omega0, T cSmago) {
        return (T)0.5 * util::sqr(cSmago*omega0*Descriptor<T>::invCs2);
    }
    static T computeOmega(T omega0, T alpha, Array<T,Descriptor<T>::d> const& j1)
    {
        T j1Norm      = norm(j1);
        T linearTerm  = alpha*j1Norm;
        T squareTerm  = (T)2*alpha*alpha*j1Norm;
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
    meta::registerTwoParamDynamics<T,Descriptor,SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor> >("SmagoAdvectionDiffusion_RLB");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>::SmagorinskyAdvectionDiffusionRLBdynamics (T omega_, T cSmago_ )
    : AdvectionDiffusionDynamics<T,Descriptor>(omega_),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>* SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>::clone() const {
    return new SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>::getId() const {
    return id;
}
 
template<typename T, template<typename U> class Descriptor>
void SmagorinskyAdvectionDiffusionRLBdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics ) 
{
    T rhoBar;
    Array<T,Descriptor<T>::d> jEq, jNeq;
    advectionDiffusionMomentTemplates<T,Descriptor>::get_rhoBar_jEq_jNeq(cell, rhoBar, jEq, jNeq);
    
    T omega0 = this->getOmega();
    T alpha = AD_SmagoOperations<T,Descriptor>::computePrefactor(omega0, cSmago);
    T omega = AD_SmagoOperations<T,Descriptor>::computeOmega(omega0, alpha, jNeq);

    T uSqr = advectionDiffusionDynamicsTemplates<T,Descriptor>::
            no_corr_rlb_collision(cell, rhoBar, jEq, jNeq, omega );
    
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

/** \param j The parameter j is defined as j = j_advDiff = rho_advDiff*u_fluid
 */
template<typename T, template<typename U> class Descriptor>
T AdvectionDiffusionBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{ 
    return advectionDiffusionDynamicsTemplates<T,Descriptor>::
            bgk_ma1_equilibrium(iPop, rhoBar, j);

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
