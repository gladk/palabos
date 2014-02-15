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

/* Orestis Malaspinas contributed this code. */

/** \file
 * MRT dynamics -- generic implementation.
 */
#ifndef EXTERNAL_FORCE_MRT_DYNAMICS_HH
#define EXTERNAL_FORCE_MRT_DYNAMICS_HH

#include "externalForceMrtDynamics.h"
#include "latticeBoltzmann/mrtTemplates.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "core/latticeStatistics.h"
#include <algorithm>
#include <limits>

namespace plb {
    
/* *************** Class GuoExternalForceMRTdynamics *********************************************** */
template<typename T, template<typename U> class Descriptor>
int GuoExternalForceMRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,GuoExternalForceMRTdynamics<T,Descriptor> >("MRT_ExternalForce_Guo");

template<typename T, template<typename U> class Descriptor>
GuoExternalForceMRTdynamics<T,Descriptor>::GuoExternalForceMRTdynamics(plint externalParam_)
    : ExternalForceDynamics<T,Descriptor>(mrtParam<T,Descriptor>().get(externalParam_).getOmega()),
      param(0),
      externalParam(externalParam_)
{ } 

template<typename T, template<typename U> class Descriptor>
GuoExternalForceMRTdynamics<T,Descriptor>::GuoExternalForceMRTdynamics(MRTparam<T,Descriptor>* param_)
    : ExternalForceDynamics<T,Descriptor>(param_->getOmega()),
      param(param_),
      externalParam(-1)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceMRTdynamics<T,Descriptor>::GuoExternalForceMRTdynamics(HierarchicUnserializer& unserializer) 
    : ExternalForceDynamics<T,Descriptor>(T()),
      param(0),
      externalParam(-1)
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceMRTdynamics<T,Descriptor>::GuoExternalForceMRTdynamics(GuoExternalForceMRTdynamics<T,Descriptor> const& rhs)
    : ExternalForceDynamics<T,Descriptor>(rhs),
      param(rhs.param ? rhs.param->clone() : 0),
      externalParam(rhs.externalParam)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceMRTdynamics<T,Descriptor>&
    GuoExternalForceMRTdynamics<T,Descriptor>::operator=(GuoExternalForceMRTdynamics<T,Descriptor> const& rhs)
{
    GuoExternalForceMRTdynamics<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceMRTdynamics<T,Descriptor>::~GuoExternalForceMRTdynamics()
{
    delete param;
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceMRTdynamics<T,Descriptor>::swap(GuoExternalForceMRTdynamics<T,Descriptor>& rhs)
{
    std::swap(param, rhs.param);
    std::swap(externalParam, rhs.externalParam);
    T tmpOmega = this->getOmega();
    this->setOmega(rhs.getOmega());
    rhs.setOmega(tmpOmega);
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceMRTdynamics<T,Descriptor>* GuoExternalForceMRTdynamics<T,Descriptor>::clone() const {
    return new GuoExternalForceMRTdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int GuoExternalForceMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceMRTdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    int useExternalParamFlag = param ? 0:1;
    serializer.addValue(useExternalParamFlag);
    if (param) {
        param->serialize(serializer);
    }
    else {
        serializer.addValue(externalParam);
    }
    ExternalForceDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceMRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
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
    ExternalForceDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceMRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    T rhoBar = momentTemplates<T,Descriptor>::get_rhoBar(cell);
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    
    Array<T,Descriptor<T>::d> u, j;
    this->computeVelocity(cell, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        j[iD] = rho * u[iD];
    }
    
    T jSqr = mrtTemp::mrtCollisionWithForce(cell, rhoBar, j, parameter->getInvM(), (T)1);
    
    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr*invRho*invRho);
//         gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T GuoExternalForceMRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor> const& GuoExternalForceMRTdynamics<T,Descriptor>::getMrtParameter() const {
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    return *parameter;
}

/* *************** Class VariableOmegaForcedMRTdynamics *********************************************** */

// Set omega to 1, because the actual multiplication by omega will be done
// during collision.
template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor> VariableOmegaForcedMRTdynamics<T,Descriptor>::param(1.);

template<typename T, template<typename U> class Descriptor>
int VariableOmegaForcedMRTdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,VariableOmegaForcedMRTdynamics<T,Descriptor> >("VariableOmegaMRT");

template<typename T, template<typename U> class Descriptor>
VariableOmegaForcedMRTdynamics<T,Descriptor>::VariableOmegaForcedMRTdynamics(T omega)
    : ExternalForceDynamics<T,Descriptor>(omega)
{ } 

template<typename T, template<typename U> class Descriptor>
VariableOmegaForcedMRTdynamics<T,Descriptor>* VariableOmegaForcedMRTdynamics<T,Descriptor>::clone() const {
    return new VariableOmegaForcedMRTdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int VariableOmegaForcedMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void VariableOmegaForcedMRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    T rhoBar = momentTemplates<T,Descriptor>::get_rhoBar(cell);
    
    Array<T,Descriptor<T>::d> u, j;
    this->computeVelocity(cell, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        j[iD] = rho * u[iD];
    }
    
    T jSqr = mrtTemp::variableOmegaMrtCollisionWithForce(cell, rhoBar, j, param.getInvM(), (T)1, this->getOmega());
    
    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr*invRho*invRho);
    }
}

template<typename T, template<typename U> class Descriptor>
T VariableOmegaForcedMRTdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium
               (iPop, rhoBar, invRho, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor> const& VariableOmegaForcedMRTdynamics<T,Descriptor>::getMrtParameter() const {
    return param;
}




/* *************** Class GuoExternalForceSmagorinskyMRTdynamics *********************************************** */
template<typename T, template<typename U> class Descriptor>
int GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor> >("MRT_ExternalForce_Guo_Smagorinsky");

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::GuoExternalForceSmagorinskyMRTdynamics(plint externalParam_, T cSmago_)
    : ExternalForceDynamics<T,Descriptor>(mrtParam<T,Descriptor>().get(externalParam_).getOmega()),
      param(0),
      externalParam(externalParam_),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::GuoExternalForceSmagorinskyMRTdynamics(MRTparam<T,Descriptor>* param_, T cSmago_)
    : ExternalForceDynamics<T,Descriptor>(param_->getOmega()),
      param(param_),
      externalParam(-1),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::GuoExternalForceSmagorinskyMRTdynamics(HierarchicUnserializer& unserializer)
    : ExternalForceDynamics<T,Descriptor>(T()),
      param(0),
      externalParam(-1)
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::GuoExternalForceSmagorinskyMRTdynamics(GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor> const& rhs)
    : ExternalForceDynamics<T,Descriptor>(rhs),
      param(rhs.param ? rhs.param->clone() : 0),
      externalParam(rhs.externalParam),
      cSmago(rhs.cSmago)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>&
    GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::operator=(GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor> const& rhs)
{
    GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::~GuoExternalForceSmagorinskyMRTdynamics()
{
    delete param;
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::swap(GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>& rhs)
{
    std::swap(param, rhs.param);
    std::swap(externalParam, rhs.externalParam);
    T tmpOmega = this->getOmega();
    this->setOmega(rhs.getOmega());
    rhs.setOmega(tmpOmega);
    std::swap(cSmago,rhs.cSmago_);
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>* GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::clone() const {
    return new GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
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
    ExternalForceDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
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
    ExternalForceDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::collide (
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

    T jSqr = mrtTemp::mrtCollisionWithForce(cell, rhoBar, j, invM_S, (T)1);

    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr*invRho*invRho);
//         gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor> const& GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>::getMrtParameter() const {
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    return *parameter;
}


/* *************** Class GuoExternalForceSmagorinskyQuasiIncMRTdynamics *********************************************** */
template<typename T, template<typename U> class Descriptor>
int GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor> >("QuasiInc_MRT_ExternalForce_Guo_Smagorinsky");

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::GuoExternalForceSmagorinskyQuasiIncMRTdynamics(plint externalParam_, T cSmago_)
    : ExternalForceDynamics<T,Descriptor>(mrtParam<T,Descriptor>().get(externalParam_).getOmega()),
      param(0),
      externalParam(externalParam_),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::GuoExternalForceSmagorinskyQuasiIncMRTdynamics(MRTparam<T,Descriptor>* param_, T cSmago_)
    : ExternalForceDynamics<T,Descriptor>(param_->getOmega()),
      param(param_),
      externalParam(-1),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::GuoExternalForceSmagorinskyQuasiIncMRTdynamics(HierarchicUnserializer& unserializer)
    : ExternalForceDynamics<T,Descriptor>(T()),
      param(0),
      externalParam(-1)
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::GuoExternalForceSmagorinskyQuasiIncMRTdynamics(GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor> const& rhs)
    : ExternalForceDynamics<T,Descriptor>(rhs),
      param(rhs.param ? rhs.param->clone() : 0),
      externalParam(rhs.externalParam),
      cSmago(rhs.cSmago)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>&
    GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::operator=(GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor> const& rhs)
{
    GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::~GuoExternalForceSmagorinskyQuasiIncMRTdynamics()
{
    delete param;
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::swap(GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>& rhs)
{
    std::swap(param, rhs.param);
    std::swap(externalParam, rhs.externalParam);
    T tmpOmega = this->getOmega();
    this->setOmega(rhs.getOmega());
    rhs.setOmega(tmpOmega);
    std::swap(cSmago,rhs.cSmago_);
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>* GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::clone() const {
    return new GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
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
    ExternalForceDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
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
    ExternalForceDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::collide (
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
        Array<T,SymmetricTensor<T,Descriptor>::n> S = (-rho*tau*Descriptor<T>::cs2+sqrt(util::sqr(rho*tau*Descriptor<T>::cs2)+normPiNeq)) / normPiNeq * PiNeq;

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

    T jSqr = mrtTemp::quasiIncMrtCollisionWithForce(cell, rhoBar, j, invM_S, (T)1);
    
    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr*invRho*invRho);
//         gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor> const& GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::getMrtParameter() const {
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    return *parameter;
}

/* *************** Class GuoExternalForceConsistentSmagorinskyMRTdynamics *********************************************** */
template<typename T, template<typename U> class Descriptor>
int GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor> >("MRT_ExternalForce_Guo_Consistent_Smagorinsky");

template<typename T, template<typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::GuoExternalForceConsistentSmagorinskyMRTdynamics(plint externalParam_, T cSmago_)
    : ExternalForceDynamics<T,Descriptor>(mrtParam<T,Descriptor>().get(externalParam_).getOmega()),
      param(0),
      externalParam(externalParam_),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::GuoExternalForceConsistentSmagorinskyMRTdynamics(MRTparam<T,Descriptor>* param_, T cSmago_)
    : ExternalForceDynamics<T,Descriptor>(param_->getOmega()),
      param(param_),
      externalParam(-1),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::GuoExternalForceConsistentSmagorinskyMRTdynamics(HierarchicUnserializer& unserializer)
    : ExternalForceDynamics<T,Descriptor>(T()),
      param(0),
      externalParam(-1)
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::GuoExternalForceConsistentSmagorinskyMRTdynamics(GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor> const& rhs)
    : ExternalForceDynamics<T,Descriptor>(rhs),
      param(rhs.param ? rhs.param->clone() : 0),
      externalParam(rhs.externalParam),
      cSmago(rhs.cSmago)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>&
    GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::operator=(GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor> const& rhs)
{
    GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::~GuoExternalForceConsistentSmagorinskyMRTdynamics()
{
    delete param;
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::swap(GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>& rhs)
{
    std::swap(param, rhs.param);
    std::swap(externalParam, rhs.externalParam);
    T tmpOmega = this->getOmega();
    this->setOmega(rhs.getOmega());
    rhs.setOmega(tmpOmega);
    std::swap(cSmago,rhs.cSmago_);
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>* GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::clone() const {
    return new GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
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
    ExternalForceDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
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
    ExternalForceDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::collide (
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
    
    T jSqr = mrtTemp::smagorinskyMrtCollisionWithForce(cell, rhoBar, j, S, cSmago, parameter->getInvM(), (T)1);
    
    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr*invRho*invRho);
    }
}

template<typename T, template<typename U> class Descriptor>
T GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor> const& GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>::getMrtParameter() const {
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    return *parameter;
}


/* *************** Class GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics *********************************************** */
template<typename T, template<typename U> class Descriptor>
int GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor> >("QuasiInc_MRT_ExternalForce_Guo_Consistent_Smagorinsky");

template<typename T, template<typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics(plint externalParam_, T cSmago_)
    : ExternalForceDynamics<T,Descriptor>(mrtParam<T,Descriptor>().get(externalParam_).getOmega()),
      param(0),
      externalParam(externalParam_),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics(MRTparam<T,Descriptor>* param_, T cSmago_)
    : ExternalForceDynamics<T,Descriptor>(param_->getOmega()),
      param(param_),
      externalParam(-1),
      cSmago(cSmago_)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics(HierarchicUnserializer& unserializer)
    : ExternalForceDynamics<T,Descriptor>(T()),
      param(0),
      externalParam(-1)
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics(GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor> const& rhs)
    : ExternalForceDynamics<T,Descriptor>(rhs),
      param(rhs.param ? rhs.param->clone() : 0),
      externalParam(rhs.externalParam),
      cSmago(rhs.cSmago)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>&
    GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::operator=(GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor> const& rhs)
{
    GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::~GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics()
{
    delete param;
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::swap(GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>& rhs)
{
    std::swap(param, rhs.param);
    std::swap(externalParam, rhs.externalParam);
    T tmpOmega = this->getOmega();
    this->setOmega(rhs.getOmega());
    rhs.setOmega(tmpOmega);
    std::swap(cSmago,rhs.cSmago_);
}

template<typename T, template<typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>* GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::clone() const {
    return new GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
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
    ExternalForceDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
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
    ExternalForceDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::collide (
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
    
    T jSqr = mrtTemp::quasiIncSmagorinskyMrtCollisionWithForce(cell, rhoBar, j, S, cSmago, parameter->getInvM(), (T)1);
    
    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr*invRho*invRho);
    }
}

template<typename T, template<typename U> class Descriptor>
T GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, (T)1, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor> const& GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>::getMrtParameter() const {
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    return *parameter;
}

}
#endif

