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
#ifndef MRT_DYNAMICS_HH
#define MRT_DYNAMICS_HH

#include "latticeBoltzmann/mrtTemplates.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "core/latticeStatistics.h"
#include <algorithm>
#include <limits>

namespace plb {

/* *************** Class MRTparam *********************************************** */

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor>::MRTparam()
{
    iniRelaxationVector();
    precompute_invM_S();
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor>::MRTparam(T omega_)
{
    iniRelaxationVector();
    for (plint iPop  = 0; iPop < Descriptor<T>::shearIndexes; ++iPop) {
        s[Descriptor<T>::shearViscIndexes[iPop]] = omega_;
    }
    precompute_invM_S();
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor>::MRTparam(T omega_, T lambda_)
{
    iniRelaxationVector();
    for (plint iPop  = 0; iPop < Descriptor<T>::shearIndexes; ++iPop) {
        s[Descriptor<T>::shearViscIndexes[iPop]] = omega_;
    }
    s[Descriptor<T>::bulkViscIndex] = lambda_;
    precompute_invM_S();
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor>::MRTparam(std::vector<T> s_)
    : s(s_)
{
    PLB_ASSERT(s.size() == Descriptor<T>::q);
    precompute_invM_S();
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor>::MRTparam(HierarchicUnserializer& unserializer)
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void MRTparam<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue((int)s.size());
    serializer.addValues(s);
}
    
template<typename T, template<typename U> class Descriptor>
void MRTparam<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    s.resize(unserializer.readValue<int>());
    unserializer.readValues(s);
    precompute_invM_S();
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor>* MRTparam<T,Descriptor>::clone() const {
    return new MRTparam<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
T MRTparam<T,Descriptor>::getParameter(plint whichParameter) const 
{
    switch (whichParameter) {
        case dynamicParams::omega_shear   : return getOmega();
        case dynamicParams::omega_bulk    : return getLambda();
        case dynamicParams::omega_epsilon : return getOmegaEpsilon();
        case dynamicParams::omega_q       : return getOmegaQ();
    };
    return 0.;
}

template<typename T, template<typename U> class Descriptor>
void MRTparam<T,Descriptor>::setParameter(plint whichParameter, T value)
{
    switch (whichParameter) {
        case dynamicParams::omega_shear   : setOmega(value);
        case dynamicParams::omega_bulk    : setLambda(value);
        case dynamicParams::omega_q       : setOmegaQ(value);
        case dynamicParams::omega_epsilon : setOmegaEpsilon(value);
    };
    precompute_invM_S();
}


template<typename T, template<typename U> class Descriptor>
std::vector<T> const& MRTparam<T,Descriptor>::getS() const 
{
    return s;
}

template<typename T, template<typename U> class Descriptor>
typename MRTparam<T,Descriptor>::MatT& MRTparam<T,Descriptor>::getInvM()
{
    return invM_S;
}

template<typename T, template<typename U> class Descriptor>
T MRTparam<T,Descriptor>::getOmega() const 
{
    return s[Descriptor<T>::shearViscIndexes[0]];
}

template<typename T, template<typename U> class Descriptor>
void MRTparam<T,Descriptor>::setOmega(T omega_) 
{
    for (plint iPop  = 0; iPop < Descriptor<T>::shearIndexes; ++iPop) {
        s[Descriptor<T>::shearViscIndexes[iPop]] = omega_;
    }
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        for (plint jPop_ = 0; jPop_ < Descriptor<T>::shearIndexes; ++jPop_) {
            plint jPop = Descriptor<T>::shearViscIndexes[jPop_];
            invM_S[iPop][jPop] = Descriptor<T>::invM[iPop][jPop] * s[jPop];
        }
    }
}

template<typename T, template<typename U> class Descriptor>
T MRTparam<T,Descriptor>::getOmegaQ() const 
{
    return s[Descriptor<T>::qViscIndexes[0]];
}

template<typename T, template<typename U> class Descriptor>
void MRTparam<T,Descriptor>::setOmegaQ(T q_) 
{
    for (int iPop = 0; iPop < Descriptor<T>::qIndexes; ++iPop) {
        s[Descriptor<T>::qViscIndexes[iPop]] = q_;
    }
    precompute_invM_S();
}

template<typename T, template<typename U> class Descriptor>
T MRTparam<T,Descriptor>::getOmegaEpsilon() const 
{
    return s[Descriptor<T>::epsilonIndex];
}

template<typename T, template<typename U> class Descriptor>
void MRTparam<T,Descriptor>::setOmegaEpsilon(T epsilon_) 
{
    s[Descriptor<T>::epsilonIndex] = epsilon_;
    precompute_invM_S();
}

template<typename T, template<typename U> class Descriptor>
T MRTparam<T,Descriptor>::getLambda() const 
{
    return s[Descriptor<T>::bulkViscIndex];
}

template<typename T, template<typename U> class Descriptor>
void MRTparam<T,Descriptor>::setLambda(T lambda_) 
{
    s[Descriptor<T>::bulkViscIndex] = lambda_;
    precompute_invM_S();
}

template<typename T, template<typename U> class Descriptor>
void MRTparam<T,Descriptor>::iniRelaxationVector()
{
    s.resize(Descriptor<T>::q);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        s[iPop] = Descriptor<T>::S[iPop];
    }
}

template<typename T, template<typename U> class Descriptor>
void MRTparam<T,Descriptor>::precompute_invM_S()
{
    PLB_ASSERT(s.size() == Descriptor<T>::q);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        for (plint jPop = 0; jPop < Descriptor<T>::q; ++jPop) {
            invM_S[iPop][jPop] = Descriptor<T>::invM[iPop][jPop] * s[jPop];
        }
    }
}


/* *************** Class MRTparamList ********************************************** */

template<typename T, template<typename U> class Descriptor>
void MRTparamList<T,Descriptor>::set(plint id, MRTparam<T,Descriptor> const& param)
{
    PLB_PRECONDITION(id>=0);
    parameters.resize(id+1);
    parameters[id] = param;
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor>&  MRTparamList<T,Descriptor>::get(plint id) {
    PLB_PRECONDITION(id>=0 && id<(plint)parameters.size());
    return parameters[id];
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor> const&  MRTparamList<T,Descriptor>::get(plint id) const
{
    PLB_PRECONDITION(id>=0 && id<(plint)parameters.size());
    return parameters[id];
}

template<typename T, template<typename U> class Descriptor>
MRTparamList<T,Descriptor>& mrtParam() {
    static MRTparamList<T,Descriptor> paramList;
    return paramList;
}


/* *************** Class MRTdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
int MRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,MRTdynamics<T,Descriptor> >("MRT");

template<typename T, template<typename U> class Descriptor>
MRTdynamics<T,Descriptor>::MRTdynamics(plint externalParam_)
    : IsoThermalBulkDynamics<T,Descriptor>(mrtParam<T,Descriptor>().get(externalParam_).getOmega()),
      param(0),
      externalParam(externalParam_)
{ } 

template<typename T, template<typename U> class Descriptor>
MRTdynamics<T,Descriptor>::MRTdynamics(MRTparam<T,Descriptor>* param_)
    : IsoThermalBulkDynamics<T,Descriptor>(param_->getOmega()),
      param(param_),
      externalParam(-1)
{ }

template<typename T, template<typename U> class Descriptor>
MRTdynamics<T,Descriptor>::MRTdynamics(HierarchicUnserializer& unserializer) 
    : IsoThermalBulkDynamics<T,Descriptor>(T()),
      param(0),
      externalParam(-1)
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
MRTdynamics<T,Descriptor>::MRTdynamics(MRTdynamics<T,Descriptor> const& rhs)
    : IsoThermalBulkDynamics<T,Descriptor>(rhs),
      param(rhs.param ? rhs.param->clone() : 0),
      externalParam(rhs.externalParam)
{ }

template<typename T, template<typename U> class Descriptor>
MRTdynamics<T,Descriptor>&
    MRTdynamics<T,Descriptor>::operator=(MRTdynamics<T,Descriptor> const& rhs)
{
    MRTdynamics<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
MRTdynamics<T,Descriptor>::~MRTdynamics()
{
    delete param;
}

template<typename T, template<typename U> class Descriptor>
void MRTdynamics<T,Descriptor>::setOmega(T omega_) {
    if (param) {
        param->setOmega(omega_);
    }
}

template<typename T, template<typename U> class Descriptor>
void MRTdynamics<T,Descriptor>::swap(MRTdynamics<T,Descriptor>& rhs)
{
    std::swap(param, rhs.param);
    std::swap(externalParam, rhs.externalParam);
    T tmpOmega = this->getOmega();
    this->setOmega(rhs.getOmega());
    rhs.setOmega(tmpOmega);
}

template<typename T, template<typename U> class Descriptor>
MRTdynamics<T,Descriptor>* MRTdynamics<T,Descriptor>::clone() const {
    return new MRTdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int MRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void MRTdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    int useExternalParamFlag = param ? 0:1;
    serializer.addValue(useExternalParamFlag);
    if (param) {
        param->serialize(serializer);
    }
    else {
        serializer.addValue(externalParam);
    }
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void MRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
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
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void MRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar,j);
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    T jSqr = mrtTemp::mrtCollision(cell, rhoBar, j, parameter->getInvM());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar) );
    }
}

template<typename T, template<typename U> class Descriptor>
void MRTdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T thetaBar, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    T jSqr = mrtTemp::mrtCollision(cell, rhoBar, j, parameter->getInvM());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar) );
    }
}

template<typename T, template<typename U> class Descriptor>
T MRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor> const& MRTdynamics<T,Descriptor>::getMrtParameter() const {
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    return *parameter;
}


/* *************** Class VariableOmegaMRTdynamics *********************************************** */

// Set omega to 1, because the actual multiplication by omega will be done
// during collision.
template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor> VariableOmegaMRTdynamics<T,Descriptor>::param(1.);

template<typename T, template<typename U> class Descriptor>
int VariableOmegaMRTdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,VariableOmegaMRTdynamics<T,Descriptor> >("VariableOmegaMRT");

template<typename T, template<typename U> class Descriptor>
VariableOmegaMRTdynamics<T,Descriptor>::VariableOmegaMRTdynamics(T omega)
    : IsoThermalBulkDynamics<T,Descriptor>(omega)
{ } 

template<typename T, template<typename U> class Descriptor>
VariableOmegaMRTdynamics<T,Descriptor>* VariableOmegaMRTdynamics<T,Descriptor>::clone() const {
    return new VariableOmegaMRTdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int VariableOmegaMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void VariableOmegaMRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar,j);
    T jSqr = mrtTemp::variableOmegaMrtCollision(cell, rhoBar, j, param.getInvM(), this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar) );
    }
}

template<typename T, template<typename U> class Descriptor>
void VariableOmegaMRTdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T thetaBar, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    T jSqr = mrtTemp::variableOmegaMrtCollision(cell, rhoBar, j, param.getInvM(), this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar) );
    }
}

template<typename T, template<typename U> class Descriptor>
T VariableOmegaMRTdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium
               (iPop, rhoBar, invRho, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
MRTparam<T,Descriptor> const& VariableOmegaMRTdynamics<T,Descriptor>::getMrtParameter() const {
    return param;
}


/* *************** Class ExternalVelocityMRTdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
int ExternalVelocityMRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,ExternalVelocityMRTdynamics<T,Descriptor> >("ExternalVelocityMRT");

template<typename T, template<typename U> class Descriptor>
ExternalVelocityMRTdynamics<T,Descriptor>::ExternalVelocityMRTdynamics(plint externalParam_)
    : IsoThermalBulkDynamics<T,Descriptor>(mrtParam<T,Descriptor>().get(externalParam_).getOmega()),
      param(0),
      externalParam(externalParam_)
{ } 

template<typename T, template<typename U> class Descriptor>
ExternalVelocityMRTdynamics<T,Descriptor>::ExternalVelocityMRTdynamics(MRTparam<T,Descriptor>* param_)
    : IsoThermalBulkDynamics<T,Descriptor>(param_->getOmega()),
      param(param_),
      externalParam(-1)
{ }

template<typename T, template<typename U> class Descriptor>
ExternalVelocityMRTdynamics<T,Descriptor>::ExternalVelocityMRTdynamics(HierarchicUnserializer& unserializer) 
    : IsoThermalBulkDynamics<T,Descriptor>(T()),
      param(0),
      externalParam(-1)
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
ExternalVelocityMRTdynamics<T,Descriptor>::ExternalVelocityMRTdynamics(ExternalVelocityMRTdynamics<T,Descriptor> const& rhs)
    : IsoThermalBulkDynamics<T,Descriptor>(rhs),
      param(rhs.param ? rhs.param->clone() : 0),
      externalParam(rhs.externalParam)
{ }

template<typename T, template<typename U> class Descriptor>
ExternalVelocityMRTdynamics<T,Descriptor>&
    ExternalVelocityMRTdynamics<T,Descriptor>::operator=(ExternalVelocityMRTdynamics<T,Descriptor> const& rhs)
{
    ExternalVelocityMRTdynamics<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
ExternalVelocityMRTdynamics<T,Descriptor>::~ExternalVelocityMRTdynamics()
{
    delete param;
}

template<typename T, template<typename U> class Descriptor>
void ExternalVelocityMRTdynamics<T,Descriptor>::swap(ExternalVelocityMRTdynamics<T,Descriptor>& rhs)
{
    std::swap(param, rhs.param);
    std::swap(externalParam, rhs.externalParam);
    T tmpOmega = this->getOmega();
    this->setOmega(rhs.getOmega());
    rhs.setOmega(tmpOmega);
}

template<typename T, template<typename U> class Descriptor>
ExternalVelocityMRTdynamics<T,Descriptor>* ExternalVelocityMRTdynamics<T,Descriptor>::clone() const {
    return new ExternalVelocityMRTdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int ExternalVelocityMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void ExternalVelocityMRTdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    int useExternalParamFlag = param ? 0:1;
    serializer.addValue(useExternalParamFlag);
    if (param) {
        param->serialize(serializer);
    }
    else {
        serializer.addValue(externalParam);
    }
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void ExternalVelocityMRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
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
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void ExternalVelocityMRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    T rhoBar = momentTemplates<T,Descriptor>::get_rhoBar(cell);
    Array<T,Descriptor<T>::d> j;
    j.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    j *= Descriptor<T>::fullRho(rhoBar);
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    T jSqr = mrtTemp::mrtCollision(cell, rhoBar, j, parameter->getInvM());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar) );
    }
}

template<typename T, template<typename U> class Descriptor>
T ExternalVelocityMRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
T ExternalVelocityMRTdynamics<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const {
    return momentTemplates<T,Descriptor>::get_rhoBar(cell);
}

template<typename T, template<typename U> class Descriptor>
void ExternalVelocityMRTdynamics<T,Descriptor>::computeRhoBarJ(Cell<T,Descriptor> const& cell,
                    T& rhoBar, Array<T,Descriptor<T>::d>& j) const
{
    rhoBar = momentTemplates<T,Descriptor>::get_rhoBar(cell);
    j.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    j *= Descriptor<T>::fullRho(rhoBar);
}

template<typename T, template<typename U> class Descriptor>
void ExternalVelocityMRTdynamics<T,Descriptor>::computeVelocity( Cell<T,Descriptor> const& cell, Array<T,Descriptor<T>::d>& u ) const {
    u.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
}

/* *************** Class IncMRTdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
int IncMRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,IncMRTdynamics<T,Descriptor> >("IncMRT");

template<typename T, template<typename U> class Descriptor>
IncMRTdynamics<T,Descriptor>::IncMRTdynamics(plint externalParam_)
    : IsoThermalBulkDynamics<T,Descriptor>(mrtParam<T,Descriptor>().get(externalParam_).getOmega()),
      param(0),
      externalParam(externalParam_)
{ } 

template<typename T, template<typename U> class Descriptor>
IncMRTdynamics<T,Descriptor>::IncMRTdynamics(MRTparam<T,Descriptor>* param_)
    : IsoThermalBulkDynamics<T,Descriptor>(param_->getOmega()),
      param(param_),
      externalParam(-1)
{ }

template<typename T, template<typename U> class Descriptor>
IncMRTdynamics<T,Descriptor>::IncMRTdynamics(HierarchicUnserializer& unserializer) 
    : IsoThermalBulkDynamics<T,Descriptor>(T()),
      param(0),
      externalParam(-1)
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
IncMRTdynamics<T,Descriptor>::IncMRTdynamics(IncMRTdynamics<T,Descriptor> const& rhs)
    : IsoThermalBulkDynamics<T,Descriptor>(rhs),
      param(rhs.param ? rhs.param->clone() : 0),
      externalParam(rhs.externalParam)
{ }

template<typename T, template<typename U> class Descriptor>
IncMRTdynamics<T,Descriptor>&
    IncMRTdynamics<T,Descriptor>::operator=(IncMRTdynamics<T,Descriptor> const& rhs)
{
    IncMRTdynamics<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
IncMRTdynamics<T,Descriptor>::~IncMRTdynamics()
{
    delete param;
}

template<typename T, template<typename U> class Descriptor>
void IncMRTdynamics<T,Descriptor>::swap(IncMRTdynamics<T,Descriptor>& rhs)
{
    std::swap(param, rhs.param);
    std::swap(externalParam, rhs.externalParam);
    T tmpOmega = this->getOmega();
    this->setOmega(rhs.getOmega());
    rhs.setOmega(tmpOmega);
}

template<typename T, template<typename U> class Descriptor>
IncMRTdynamics<T,Descriptor>* IncMRTdynamics<T,Descriptor>::clone() const {
    return new IncMRTdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int IncMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void IncMRTdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    int useExternalParamFlag = param ? 0:1;
    serializer.addValue(useExternalParamFlag);
    if (param) {
        param->serialize(serializer);
    }
    else {
        serializer.addValue(externalParam);
    }
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void IncMRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
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
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void IncMRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar,j);
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    T jSqr = mrtTemp::quasiIncMrtCollision(cell, rhoBar, j, parameter->getInvM());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar) );
    }
}

template<typename T, template<typename U> class Descriptor>
T IncMRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, (T)1, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
bool IncMRTdynamics<T,Descriptor>::velIsJ() const {
    return true;
}

template<typename T, template<typename U> class Descriptor>
void IncMRTdynamics<T,Descriptor>::computeVelocity( Cell<T,Descriptor> const& cell,
                                                                 Array<T,Descriptor<T>::d>& u ) const 
{
    T dummyRhoBar;
    this->computeRhoBarJ(cell, dummyRhoBar, u);
}

/* *************** Class ExternalVelocityIncMRTdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
int ExternalVelocityIncMRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,ExternalVelocityIncMRTdynamics<T,Descriptor> >("ExternalVelocityIncMRT");

template<typename T, template<typename U> class Descriptor>
ExternalVelocityIncMRTdynamics<T,Descriptor>::ExternalVelocityIncMRTdynamics(plint externalParam_)
    : IsoThermalBulkDynamics<T,Descriptor>(mrtParam<T,Descriptor>().get(externalParam_).getOmega()),
      param(0),
      externalParam(externalParam_)
{ } 

template<typename T, template<typename U> class Descriptor>
ExternalVelocityIncMRTdynamics<T,Descriptor>::ExternalVelocityIncMRTdynamics(MRTparam<T,Descriptor>* param_)
    : IsoThermalBulkDynamics<T,Descriptor>(param_->getOmega()),
      param(param_),
      externalParam(-1)
{ }

template<typename T, template<typename U> class Descriptor>
ExternalVelocityIncMRTdynamics<T,Descriptor>::ExternalVelocityIncMRTdynamics(HierarchicUnserializer& unserializer) 
    : IsoThermalBulkDynamics<T,Descriptor>(T()),
      param(0),
      externalParam(-1)
{
    unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
ExternalVelocityIncMRTdynamics<T,Descriptor>::ExternalVelocityIncMRTdynamics(
        ExternalVelocityIncMRTdynamics<T,Descriptor> const& rhs)
    : IsoThermalBulkDynamics<T,Descriptor>(rhs),
      param(rhs.param ? rhs.param->clone() : 0),
      externalParam(rhs.externalParam)
{ }

template<typename T, template<typename U> class Descriptor>
ExternalVelocityIncMRTdynamics<T,Descriptor>&
    ExternalVelocityIncMRTdynamics<T,Descriptor>::operator=(ExternalVelocityIncMRTdynamics<T,Descriptor> const& rhs)
{
    ExternalVelocityIncMRTdynamics<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
ExternalVelocityIncMRTdynamics<T,Descriptor>::~ExternalVelocityIncMRTdynamics()
{
    delete param;
}

template<typename T, template<typename U> class Descriptor>
void ExternalVelocityIncMRTdynamics<T,Descriptor>::swap(ExternalVelocityIncMRTdynamics<T,Descriptor>& rhs)
{
    std::swap(param, rhs.param);
    std::swap(externalParam, rhs.externalParam);
    T tmpOmega = this->getOmega();
    this->setOmega(rhs.getOmega());
    rhs.setOmega(tmpOmega);
}

template<typename T, template<typename U> class Descriptor>
ExternalVelocityIncMRTdynamics<T,Descriptor>* ExternalVelocityIncMRTdynamics<T,Descriptor>::clone() const {
    return new ExternalVelocityIncMRTdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int ExternalVelocityIncMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void ExternalVelocityIncMRTdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    int useExternalParamFlag = param ? 0:1;
    serializer.addValue(useExternalParamFlag);
    if (param) {
        param->serialize(serializer);
    }
    else {
        serializer.addValue(externalParam);
    }
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void ExternalVelocityIncMRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
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
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void ExternalVelocityIncMRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    T rhoBar = momentTemplates<T,Descriptor>::get_rhoBar(cell);
    Array<T,Descriptor<T>::d> j;
    j.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    MRTparam<T,Descriptor>* parameter = param ? param : &(mrtParam<T,Descriptor>().get(externalParam));
    T jSqr = mrtTemp::quasiIncMrtCollision(cell, rhoBar, j, parameter->getInvM());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar) );
    }
}

template<typename T, template<typename U> class Descriptor>
T ExternalVelocityIncMRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, (T)1, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
bool ExternalVelocityIncMRTdynamics<T,Descriptor>::velIsJ() const {
    return true;
}

template<typename T, template<typename U> class Descriptor>
T ExternalVelocityIncMRTdynamics<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const {
    return momentTemplates<T,Descriptor>::get_rhoBar(cell);
}

template<typename T, template<typename U> class Descriptor>
void ExternalVelocityIncMRTdynamics<T,Descriptor>::computeRhoBarJ(Cell<T,Descriptor> const& cell,
                    T& rhoBar, Array<T,Descriptor<T>::d>& j) const
{
    rhoBar = momentTemplates<T,Descriptor>::get_rhoBar(cell);
    j.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
}

template<typename T, template<typename U> class Descriptor>
void ExternalVelocityIncMRTdynamics<T,Descriptor>::computeVelocity( Cell<T,Descriptor> const& cell,
                                                                 Array<T,Descriptor<T>::d>& u ) const 
{
    T dummyRhoBar;
    computeRhoBarJ(cell, dummyRhoBar, u);
}

}

#endif  // MRT_DYNAMICS_HH

