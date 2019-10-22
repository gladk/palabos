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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef DYNAMICS_HH
#define DYNAMICS_HH

#include "core/dynamics.h"
#include "core/cell.h"
#include "core/dynamicsIdentifiers.h"
#include "core/hierarchicSerializer.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"
#include "core/latticeStatistics.h"
#include "multiGrid/multiGridUtil.h"
#include <algorithm>
#include <limits>

namespace plb {

/* *************** Class Dynamics ******************************************* */

/** By default, this method returns 0 ("undefined").  */
template<typename T, template<typename U> class Descriptor>
int Dynamics<T,Descriptor>::getId() const {
    return 0;
}

/** By default, this method returns false, as velocity is j/rho. **/
template<typename T, template<typename U> class Descriptor>
bool Dynamics<T,Descriptor>::velIsJ() const {
    return false;
}

/** By default, dynamics classes are non-composite.  */
template<typename T, template<typename U> class Descriptor>
bool Dynamics<T,Descriptor>::isComposite() const {
    return false;
}

/** By default, this method yields true.  */
template<typename T, template<typename U> class Descriptor>
bool Dynamics<T,Descriptor>::isComposeable() const {
    return true;
}

/** By default, this method yields false.  */
template<typename T, template<typename U> class Descriptor>
bool Dynamics<T,Descriptor>::isBoundary() const {
    return false;
}

/** By default, this method yields false.  */
template<typename T, template<typename U> class Descriptor>
bool Dynamics<T,Descriptor>::isAdvectionDiffusion() const {
    return false;
}

/** By default, this method yields false.  */
template<typename T, template<typename U> class Descriptor>
bool Dynamics<T,Descriptor>::isEntropic() const {
    return false;
}

/** By default, this method yields true.  */
template<typename T, template<typename U> class Descriptor>
bool Dynamics<T,Descriptor>::hasMoments() const {
    return true;
}

/** By default, this method yields false.  */
template<typename T, template<typename U> class Descriptor>
bool Dynamics<T,Descriptor>::isNonLocal() const {
    return false;
}

template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::setRelaxationFrequencies(Array<T, Descriptor<T>::q> const& frequencies) {
    setOmega(frequencies[0]);
}

template<typename T, template<typename U> class Descriptor>
Array<T, Descriptor<T>::q> Dynamics<T,Descriptor>::getRelaxationFrequencies() const {
    Array<T, Descriptor<T>::q> frequencies;
    for (plint i=0; i<Descriptor<T>::q; ++i) {
        frequencies[i] = getOmega();
    }
    return frequencies;
}

template<typename T, template<typename U> class Descriptor>
T Dynamics<T,Descriptor>::getParameter(plint whichParameter) const {
    if (whichParameter == dynamicParams::omega_shear) {
        return getOmega();
    }
    return 0.;
}

template<typename T, template<typename U> class Descriptor>
T Dynamics<T,Descriptor>::getDynamicParameter(plint whichParameter, Cell<T,Descriptor> const& cell) const {
    // Parameter not implemented.
    return ((T) 0);
}

template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::setParameter(plint whichParameter, T value) {
    if (whichParameter == dynamicParams::omega_shear) {
        setOmega(value);
    }
}

template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(this->getOmega());
}

template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    this->setOmega(unserializer.readValue<T>());
}

template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar,
        Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat )
{
    T oldRhoBar;
    Array<T,Descriptor<T>::d> oldJ;
    computeRhoBarJ(cell, oldRhoBar, oldJ);
    T oldJsqr = normSqr(oldJ);
    T jSqr = normSqr(j);
    for (int iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        T oldEq = computeEquilibrium(iPop, oldRhoBar, oldJ, oldJsqr);
        T newEq = computeEquilibrium(iPop, rhoBar, j, jSqr);
        cell[iPop] += newEq - oldEq;
    }
    collide(cell, stat);
}


template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::computeEquilibria (
        Array<T,Descriptor<T>::q>& fEq,  T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{
    for (int iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        fEq[iPop] = computeEquilibrium(iPop, rhoBar, j, jSqr, thetaBar);
    }
}

template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::rescale(int dxScale, int dtScale) {
    int dimDx = 2;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    T omega = this->getOmega();
    T nu_cs2 = (T)1/omega - (T)1/(T)2;
    this->setOmega( (T)1/(scaleFactor*nu_cs2+(T)1/(T)2) );
}

template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::getPopulations(Cell<T,Descriptor> const& cell, Array<T,Descriptor<T>::q>& f) const {
    f = cell.getRawPopulations();
}

template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::getExternalField (
        Cell<T,Descriptor> const& cell, plint pos, plint size, T* ext) const
{
    PLB_PRECONDITION(pos+size <= Descriptor<T>::ExternalField::numScalars);
    T const* externalData = cell.getExternal(pos);
    for (plint iExt=0; iExt<size; ++iExt) {
        ext[iExt] = externalData[iExt];
    }
}

template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::setPopulations(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::q> const& f)
{
    cell.getRawPopulations() = f;
}

template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::setExternalField (
        Cell<T,Descriptor>& cell, plint pos, plint size, const T* ext)
{
    PLB_PRECONDITION(pos+size <= Descriptor<T>::ExternalField::numScalars);
    T* externalData = cell.getExternal(pos);
    for (plint iExt=0; iExt<size; ++iExt) {
        externalData[iExt] = ext[iExt];
    }
}

/** This method does nothing by default, unless overloaded in
 *  a deriving class. This is OK, because it's the behavior
 *  users will expect: if they use "defineSomething" on a cell,
 *  they want it to define the boundary condition on boundary
 *  nodes, but nothing to happen in bulk nodes. In this way,
 *  they can be lazy and loop over the whole domain instead
 *  of tracking boundary nodes explicitly.
 **/
template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::defineDensity(Cell<T,Descriptor>& cell, T density) { }

/** This method does nothing by default, unless overloaded in
 *  a deriving class. This is OK, because it's the behavior
 *  users will expect: if they use "defineSomething" on a cell,
 *  they want it to define the boundary condition on boundary
 *  nodes, but nothing to happen in bulk nodes. In this way,
 *  they can be lazy and loop over the whole domain instead
 *  of tracking boundary nodes explicitly.
 **/
template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::defineVelocity(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& u) { }

/** This method does nothing by default, unless overloaded in
 *  a deriving class. This is OK, because it's the behavior
 *  users will expect: if they use "defineSomething" on a cell,
 *  they want it to define the boundary condition on boundary
 *  nodes, but nothing to happen in bulk nodes. In this way,
 *  they can be lazy and loop over the whole domain instead
 *  of tracking boundary nodes explicitly.
 **/
template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::defineTemperature(Cell<T,Descriptor>& cell, T temperature) { }

/** This method does nothing by default, unless overloaded in
 *  a deriving class. This is OK, because it's the behavior
 *  users will expect: if they use "defineSomething" on a cell,
 *  they want it to define the boundary condition on boundary
 *  nodes, but nothing to happen in bulk nodes. In this way,
 *  they can be lazy and loop over the whole domain instead
 *  of tracking boundary nodes explicitly.
 **/
template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::defineHeatFlux(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& q) { }

/** This method does nothing by default, unless overloaded in
 *  a deriving class. This is OK, because it's the behavior
 *  users will expect: if they use "defineSomething" on a cell,
 *  they want it to define the boundary condition on boundary
 *  nodes, but nothing to happen in bulk nodes. In this way,
 *  they can be lazy and loop over the whole domain instead
 *  of tracking boundary nodes explicitly.
 **/
template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::definePiNeq(Cell<T,Descriptor>& cell, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq)
{ }

/** This method does nothing by default, unless overloaded in
 *  a deriving class. This is OK, because it's the behavior
 *  users will expect: if they use "defineSomething" on a cell,
 *  they want it to define the boundary condition on boundary
 *  nodes, but nothing to happen in bulk nodes. In this way,
 *  they can be lazy and loop over the whole domain instead
 *  of tracking boundary nodes explicitly.
 **/
template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::defineMoment(Cell<T,Descriptor>& cell, plint momentId, T const* value)
{ }


/* *************** Class BasicBulkDynamics *************************** */

template<typename T, template<typename U> class Descriptor>
BasicBulkDynamics<T,Descriptor>::BasicBulkDynamics(T omega_)
    : omega(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
T BasicBulkDynamics<T,Descriptor>::computeDensity(Cell<T,Descriptor> const& cell) const
{
    return momentTemplates<T,Descriptor>::compute_rho(cell);
}

template<typename T, template<typename U> class Descriptor>
T BasicBulkDynamics<T,Descriptor>::computePressure(Cell<T,Descriptor> const& cell) const
{
    return Descriptor<T>::cs2 * computeDensity(cell) * computeTemperature(cell);
}

template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::computeVelocity (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& u) const
{
    momentTemplates<T,Descriptor>::compute_uLb(cell, u);
}

/** Defaults to 1.
 */
template<typename T, template<typename U> class Descriptor>
T BasicBulkDynamics<T,Descriptor>::computeTemperature(Cell<T,Descriptor> const& cell) const
{
    return (T)1;
}

/** Defaults to zero.
 */
template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::computePiNeq (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    PiNeq.resetToZero();
}

/** Defaults to zero.
 */
template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::computeShearStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& stress ) const
{
    stress.resetToZero();
}

/** Defaults to zero.
 */
template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::computeHeatFlux (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& q) const
{
    q.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::computeMoment(
        Cell<T,Descriptor> const& cell, plint momentId, T* moment ) const
{
    PLB_PRECONDITION( false );
}

template<typename T, template<typename U> class Descriptor>
T BasicBulkDynamics<T,Descriptor>::getOmega() const {
    return omega;
}

template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::setOmega(T omega_) {
    omega = omega_;
}

template<typename T, template<typename U> class Descriptor>
T BasicBulkDynamics<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const
{
    return momentTemplates<T,Descriptor>::get_rhoBar(cell);
}

template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j ) const
{
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
}

template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::computeRhoBarJPiNeq (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
}

template<typename T, template<typename U> class Descriptor>
T BasicBulkDynamics<T,Descriptor>::computeEbar(Cell<T,Descriptor> const& cell) const
{
    return momentTemplates<T,Descriptor>::get_eBar(cell);
}


/* *************** Class CompositeDynamics *************************** */

template<typename T, template<typename U> class Descriptor>
CompositeDynamics<T,Descriptor>::CompositeDynamics(Dynamics<T,Descriptor>* baseDynamics_,
                                                   bool automaticPrepareCollision_)
    : baseDynamics(baseDynamics_),
      automaticPrepareCollision(automaticPrepareCollision_)
{ }

template<typename T, template<typename U> class Descriptor>
CompositeDynamics<T,Descriptor>::CompositeDynamics(CompositeDynamics<T,Descriptor> const& rhs)
    : baseDynamics(rhs.baseDynamics->clone()),
      automaticPrepareCollision(rhs.automaticPrepareCollision)
{ }

template<typename T, template<typename U> class Descriptor>
CompositeDynamics<T,Descriptor>::~CompositeDynamics() {
    delete baseDynamics;
}

template<typename T, template<typename U> class Descriptor>
CompositeDynamics<T,Descriptor>& CompositeDynamics<T,Descriptor>::operator=(CompositeDynamics<T,Descriptor> const& rhs)
{
    delete baseDynamics;
    baseDynamics = rhs.baseDynamics->clone();
    automaticPrepareCollision = rhs.automaticPrepareCollision;
    return *this;
}

template<typename T, template<typename U> class Descriptor>
CompositeDynamics<T,Descriptor>* CompositeDynamics<T,Descriptor>::cloneWithNewBase (
        Dynamics<T,Descriptor>* baseDynamics_) const
{
    // First, create a clone, based on the dynamic type of CompositeDynamics
    CompositeDynamics<T,Descriptor>* newDynamics = clone();
    // Then, replace its original baseDynamics by the one we've received as parameter
    newDynamics->replaceBaseDynamics(baseDynamics_);
    return newDynamics;
}

template<typename T, template<typename U> class Descriptor>
bool CompositeDynamics<T,Descriptor>::isComposite() const {
    return true;
}

template<typename T, template<typename U> class Descriptor>
bool CompositeDynamics<T,Descriptor>::velIsJ() const {
    return this->getBaseDynamics().velIsJ();
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::replaceBaseDynamics(Dynamics<T,Descriptor>* newBaseDynamics) {
    delete baseDynamics;
    baseDynamics = newBaseDynamics;
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>& CompositeDynamics<T,Descriptor>::getBaseDynamics() {
    return *baseDynamics;
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor> const& CompositeDynamics<T,Descriptor>::getBaseDynamics() const {
    return *baseDynamics;
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(this->doesAutomaticPrepareCollision());
    serializer.nextDynamics(this->getBaseDynamics().getId());
    this->getBaseDynamics().serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    PLB_PRECONDITION( unserializer.getId() == this->getId() );
    this->toggleAutomaticPrepareCollision(unserializer.readValue<bool>());
    // If the composite dynamics is fully constructed, just initialize the variables
    //   from the unserializer.
    if (baseDynamics) {
        this->getBaseDynamics().unserialize(unserializer);
    }
    // If the composite dynamics is being constructed, then the base dynamics must
    //   be newly created.
    else {
        baseDynamics = meta::dynamicsRegistration<T,Descriptor>().generate(unserializer);
    }
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    if (automaticPrepareCollision) prepareCollision(cell);
    baseDynamics -> collide(cell, statistics);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar,
        Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat )
{
    if (automaticPrepareCollision) prepareCollision(cell);
    baseDynamics -> collideExternal(cell, rhoBar, j, thetaBar, stat);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar) const
{
    return baseDynamics -> computeEquilibrium(iPop, rhoBar, j, jSqr, thetaBar);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::computeEquilibria (
         Array<T,Descriptor<T>::q>& fEq, T rhoBar,
         Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{
    baseDynamics -> computeEquilibria(fEq, rhoBar,j, jSqr, thetaBar);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::regularize (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar ) const
{
    baseDynamics -> regularize(cell, rhoBar, j, jSqr, PiNeq, thetaBar);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::computeDensity(Cell<T,Descriptor> const& cell) const {
    return this->getBaseDynamics().computeDensity(cell);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::computePressure(Cell<T,Descriptor> const& cell) const {
    return this->getBaseDynamics().computePressure(cell);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::computeVelocity(Cell<T,Descriptor> const& cell,
                                                          Array<T,Descriptor<T>::d>& u ) const
{
    this->getBaseDynamics().computeVelocity(cell, u);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::computeTemperature(Cell<T,Descriptor> const& cell) const {
    return this->getBaseDynamics().computeTemperature(cell);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::computePiNeq (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    return this->getBaseDynamics().computePiNeq(cell, PiNeq);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::computeShearStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& stress ) const
{
    return this->getBaseDynamics().computeShearStress(cell, stress);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::computeHeatFlux(Cell<T,Descriptor> const& cell,
                                                          Array<T,Descriptor<T>::d>& q ) const
{
    this->getBaseDynamics().computeHeatFlux(cell, q);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::computeMoment (
        Cell<T,Descriptor> const& cell, plint momentId, T* moment ) const
{
    this->getBaseDynamics().computeMoment(cell, momentId, moment);
}

template<typename T, template<typename U> class Descriptor>
plint CompositeDynamics<T,Descriptor>::numDecomposedVariables(plint order) const {
    return baseDynamics->numDecomposedVariables(order);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::decompose (
        Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order ) const
{
    baseDynamics->decompose(cell, rawData, order);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::recompose (
        Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order ) const
{
    baseDynamics->recompose(cell, rawData, order);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::rescale (
        std::vector<T>& rawData, T xDxInv, T xDt, plint order ) const
{
    baseDynamics->rescale(rawData, xDxInv, xDt, order);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const
{
    return this->getBaseDynamics().computeRhoBar(cell);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j ) const
{
    this->getBaseDynamics().computeRhoBarJ(cell, rhoBar, j);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::computeRhoBarJPiNeq (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    this->getBaseDynamics().computeRhoBarJPiNeq(cell, rhoBar, j, PiNeq);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::computeEbar(Cell<T,Descriptor> const& cell) const {
    return this->getBaseDynamics().computeEbar(cell);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::getOmega() const {
    return baseDynamics -> getOmega();
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::setOmega(T omega_) {
    baseDynamics -> setOmega(omega_);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::getParameter(plint whichParameter) const {
    return baseDynamics -> getParameter(whichParameter);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::setParameter(plint whichParameter, T value) {
    baseDynamics -> setParameter(whichParameter, value);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::getPopulations(Cell<T,Descriptor> const& cell, Array<T,Descriptor<T>::q>& f) const {
    baseDynamics -> getPopulations(cell, f);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::getExternalField(Cell<T,Descriptor> const& cell,
                                                    plint pos, plint size, T* ext) const
{
    baseDynamics -> getExternalField(cell, pos, size, ext);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::setPopulations(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::q> const& f)
{
    baseDynamics -> setPopulations(cell, f);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::setExternalField(Cell<T,Descriptor>& cell, plint pos,
                                                    plint size, const T* ext)
{
    baseDynamics -> setExternalField(cell, pos, size, ext);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::defineDensity(Cell<T,Descriptor>& cell, T density) {
    baseDynamics -> defineDensity(cell, density);
}


template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::defineVelocity(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& u) {
    baseDynamics -> defineVelocity(cell, u);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::defineTemperature(Cell<T,Descriptor>& cell, T temperature) {
    baseDynamics -> defineTemperature(cell, temperature);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::defineHeatFlux(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& q) {
    baseDynamics -> defineHeatFlux(cell, q);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::definePiNeq(Cell<T,Descriptor>& cell,
                            Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq)
{
    baseDynamics -> definePiNeq(cell, PiNeq);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::defineMoment(Cell<T,Descriptor>& cell, plint momentId, T const* value) {
    baseDynamics -> defineMoment(cell, momentId, value);
}

template<typename T, template<typename U> class Descriptor>
bool CompositeDynamics<T,Descriptor>::doesAutomaticPrepareCollision() const {
    return automaticPrepareCollision;
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::toggleAutomaticPrepareCollision(bool flag) {
    automaticPrepareCollision = flag;
}

/* *************** Class PreparePopulationsDynamics *********************** */

template<typename T, template<typename U> class Descriptor>
PreparePopulationsDynamics<T,Descriptor>::PreparePopulationsDynamics (
        Dynamics<T,Descriptor>* baseDynamics_,
        bool automaticPrepareCollision_)
    : CompositeDynamics<T,Descriptor>(baseDynamics_, automaticPrepareCollision_)
{ }

template<typename T, template<typename U> class Descriptor>
void PreparePopulationsDynamics<T,Descriptor>::prepareCollision(Cell<T,Descriptor>& cell) {
    completePopulations(cell);
}


/* *************** Class BulkCompositeDynamics *********************** */

template<typename T, template<typename U> class Descriptor>
BulkCompositeDynamics<T,Descriptor>::BulkCompositeDynamics (
        Dynamics<T,Descriptor>* baseDynamics_, bool automaticPrepareCollision_)
    : PreparePopulationsDynamics<T,Descriptor>(baseDynamics_, automaticPrepareCollision_)
{ }



/* *************** Class BounceBack ********************************** */

template<typename T, template<typename U> class Descriptor>
int BounceBack<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,BounceBack<T,Descriptor> >("BounceBack");

template<typename T, template<typename U> class Descriptor>
BounceBack<T,Descriptor>::BounceBack(T rho_)
    : rho(rho_)
{ }

template<typename T, template<typename U> class Descriptor>
BounceBack<T,Descriptor>::BounceBack(HierarchicUnserializer& unserializer)
    : rho(T())
{
    this->unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
BounceBack<T,Descriptor>* BounceBack<T,Descriptor>::clone() const {
    return new BounceBack<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int BounceBack<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    Dynamics<T,Descriptor>::serialize(serializer);
    serializer.addValue(rho);
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer) 
{
    PLB_PRECONDITION( unserializer.getId() == this->getId() );
    Dynamics<T,Descriptor>::unserialize(unserializer);
    rho = unserializer.readValue<T>();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    for (plint iPop=1; iPop <= Descriptor<T>::q/2; ++iPop) {
        std::swap(cell[iPop], cell[iPop+Descriptor<T>::q/2]);
    }
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar,
        Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat )
{
    collide(cell, stat);
}

template<typename T, template<typename U> class Descriptor>
T BounceBack<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                               T jSqr, T thetaBar) const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::regularize(
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar ) const
{ }

template<typename T, template<typename U> class Descriptor>
T BounceBack<T,Descriptor>::computeDensity(Cell<T,Descriptor> const& cell) const {
    return rho;
}

template<typename T, template<typename U> class Descriptor>
T BounceBack<T,Descriptor>::computePressure(Cell<T,Descriptor> const& cell) const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::computeVelocity (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& u) const
{
    u.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
T BounceBack<T,Descriptor>::computeTemperature(Cell<T,Descriptor> const& cell) const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::computePiNeq (
        Cell<T,Descriptor> const& cell,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    PiNeq.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::computeShearStress (
        Cell<T,Descriptor> const& cell,
        Array<T,SymmetricTensor<T,Descriptor>::n>& stress ) const
{
    stress.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::computeHeatFlux (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& q) const
{
    q.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::computeMoment (
        Cell<T,Descriptor> const& cell, plint momentId, T* moment) const
{ }

template<typename T, template<typename U> class Descriptor>
T BounceBack<T,Descriptor>::getOmega() const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::setOmega(T omega_)
{ }

template<typename T, template<typename U> class Descriptor>
T BounceBack<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const {
    return Descriptor<T>::rhoBar(rho);
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j) const
{
    rhoBar = Descriptor<T>::rhoBar(rho);
    j.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::computeRhoBarJPiNeq (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    rhoBar = Descriptor<T>::rhoBar(rho);
    j.resetToZero();
    PiNeq.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
T BounceBack<T,Descriptor>::computeEbar(Cell<T,Descriptor> const& cell) const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
plint BounceBack<T,Descriptor>::numDecomposedVariables(plint order) const {
    return Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars;
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::decompose (
        Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order ) const
{
    rawData.resize(numDecomposedVariables(order));
    std::fill(rawData.begin(), rawData.end(), T());
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::recompose (
        Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order ) const
{
//     PLB_PRECONDITION( (plint)rawData.size() == numDecomposedVariables(order) );
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::rescale (
        std::vector<T>& rawData, T xDxInv, T xDt, plint order ) const
{ }

template<typename T, template<typename U> class Descriptor>
bool BounceBack<T,Descriptor>::isBoundary() const {
    return true;
}

template<typename T, template<typename U> class Descriptor>
bool BounceBack<T,Descriptor>::hasMoments() const {
    return false;
}


/* *************** Class SpecularReflection ********************************** */

template<typename T, template<typename U> class Descriptor>
int SpecularReflection<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,SpecularReflection<T,Descriptor> >("SpecularReflection");

template<typename T, template<typename U> class Descriptor>
SpecularReflection<T,Descriptor>::SpecularReflection(T rho_)
    : rho(rho_)
{
    // Configure the reflection vector so that SpecularReflection behaves like BounceBack.
    for (int iD = 0; iD < Descriptor<T>::d; iD++) {
        reflection[iD] = -1;
    }
}

template<typename T, template<typename U> class Descriptor>
SpecularReflection<T,Descriptor>::SpecularReflection(Array<bool,Descriptor<T>::d> const& reflectOnPlaneNormalToAxis, T rho_)
    : rho(rho_)
{
#ifdef PLB_DEBUG
    bool reflectSomewhere = false;
    for (int iD = 0; iD < Descriptor<T>::d; iD++) {
        reflectSomewhere = reflectSomewhere || reflectOnPlaneNormalToAxis[iD];
    }
#endif
    PLB_ASSERT(reflectSomewhere);

    Array<int,Descriptor<T>::d> reflectNormalToAxis;
    for (int iD = 0; iD < Descriptor<T>::d; iD++) {
        reflectNormalToAxis[iD] = reflectOnPlaneNormalToAxis[iD] ? 1 : 0;
    }

    reflection.resetToZero();
    for (int iD = 0; iD < Descriptor<T>::d; iD++) {
        int reflect = 0;
        for (int jD = 0; jD < Descriptor<T>::d; jD++) {
            if (iD != jD) {
                reflect = reflect || reflectNormalToAxis[jD];
            }
        }
        reflection[iD] = reflect || reflection[iD];
    }

    for (int iD = 0; iD < Descriptor<T>::d; iD++) {
        reflection[iD] = reflection[iD] == 1 ? -1 : 1;
    }
}

template<typename T, template<typename U> class Descriptor>
SpecularReflection<T,Descriptor>::SpecularReflection(HierarchicUnserializer& unserializer)
    : rho(T())
{
    // Configure the reflection vector so that SpecularReflection behaves like BounceBack.
    for (int iD = 0; iD < Descriptor<T>::d; iD++) {
        reflection[iD] = -1;
    }
    this->unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
SpecularReflection<T,Descriptor>* SpecularReflection<T,Descriptor>::clone() const
{
    return new SpecularReflection<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int SpecularReflection<T,Descriptor>::getId() const
{
    return id;
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    Dynamics<T,Descriptor>::serialize(serializer);
    serializer.addValue(rho);

    // We transform the "reflection" vector into the "reflectionInfo" vector,
    // to minimize the data serialized (from int to char).
    Array<char,Descriptor<T>::d> reflectionInfo;
    for (int iD = 0; iD < Descriptor<T>::d; iD++) {
        reflectionInfo[iD] = reflection[iD] == 1 ? 1 : 0;   // The "reflection" vector can only have the values 1, and -1.
    }
    serializer.addValues<char,Descriptor<T>::d>(reflectionInfo);
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer) 
{
    PLB_PRECONDITION( unserializer.getId() == this->getId() );
    Dynamics<T,Descriptor>::unserialize(unserializer);
    rho = unserializer.readValue<T>();

    Array<char,Descriptor<T>::d> reflectionInfo;
    unserializer.readValues<char,Descriptor<T>::d>(reflectionInfo);
    for (int iD = 0; iD < Descriptor<T>::d; iD++) {
        reflection[iD] = reflectionInfo[iD] == 1 ? 1 : -1;  // The "reflectionInfo" vector can only have the values 0 and 1.
    }
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    Array<int,Descriptor<T>::q> swapped;
    swapped.resetToZero();
    for (plint iPop = 1; iPop < Descriptor<T>::q; iPop++) {
        if (!swapped[iPop]) {
            Array<int,Descriptor<T>::d> v; 
            for (plint iD = 0; iD < Descriptor<T>::d; iD++) {
                v[iD] = reflection[iD] * Descriptor<T>::c[indexTemplates::opposite<Descriptor<T> >(iPop)][iD];
            }

            plint jPop = indexTemplates::findVelocity<Descriptor<T> >(v);
            PLB_ASSERT(jPop != Descriptor<T>::q);

            std::swap(cell[iPop], cell[jPop]);
            swapped[iPop] = 1;
            swapped[jPop] = 1;
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar,
        Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat )
{
    collide(cell, stat);
}

template<typename T, template<typename U> class Descriptor>
T SpecularReflection<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                               T jSqr, T thetaBar) const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::regularize(
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar ) const
{ }

template<typename T, template<typename U> class Descriptor>
T SpecularReflection<T,Descriptor>::computeDensity(Cell<T,Descriptor> const& cell) const
{
    return rho;
}

template<typename T, template<typename U> class Descriptor>
T SpecularReflection<T,Descriptor>::computePressure(Cell<T,Descriptor> const& cell) const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::computeVelocity (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& u) const
{
    u.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
T SpecularReflection<T,Descriptor>::computeTemperature(Cell<T,Descriptor> const& cell) const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::computePiNeq (
        Cell<T,Descriptor> const& cell,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    PiNeq.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::computeShearStress (
        Cell<T,Descriptor> const& cell,
        Array<T,SymmetricTensor<T,Descriptor>::n>& stress ) const
{
    stress.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::computeHeatFlux (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& q) const
{
    q.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::computeMoment (
        Cell<T,Descriptor> const& cell, plint momentId, T* moment) const
{ }

template<typename T, template<typename U> class Descriptor>
T SpecularReflection<T,Descriptor>::getOmega() const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::setOmega(T omega_)
{ }

template<typename T, template<typename U> class Descriptor>
T SpecularReflection<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const
{
    return Descriptor<T>::rhoBar(rho);
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j) const
{
    rhoBar = Descriptor<T>::rhoBar(rho);
    j.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::computeRhoBarJPiNeq (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    rhoBar = Descriptor<T>::rhoBar(rho);
    j.resetToZero();
    PiNeq.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
T SpecularReflection<T,Descriptor>::computeEbar(Cell<T,Descriptor> const& cell) const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
plint SpecularReflection<T,Descriptor>::numDecomposedVariables(plint order) const
{
    return Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars;
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::decompose (
        Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order ) const
{
    rawData.resize(numDecomposedVariables(order));
    std::fill(rawData.begin(), rawData.end(), T());
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::recompose (
        Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order ) const
{
//     PLB_PRECONDITION( (plint)rawData.size() == numDecomposedVariables(order) );
}

template<typename T, template<typename U> class Descriptor>
void SpecularReflection<T,Descriptor>::rescale (
        std::vector<T>& rawData, T xDxInv, T xDt, plint order ) const
{ }

template<typename T, template<typename U> class Descriptor>
bool SpecularReflection<T,Descriptor>::isBoundary() const
{
    return true;
}

template<typename T, template<typename U> class Descriptor>
bool SpecularReflection<T,Descriptor>::hasMoments() const
{
    return false;
}


/* *************** Class NoDynamics ********************************** */

template<typename T, template<typename U> class Descriptor>
int NoDynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,NoDynamics<T,Descriptor> >("NoDynamics");

template<typename T, template<typename U> class Descriptor>
NoDynamics<T,Descriptor>::NoDynamics(T rho_)
    : rho(rho_)
{ }

template<typename T, template<typename U> class Descriptor>
NoDynamics<T,Descriptor>::NoDynamics(HierarchicUnserializer& unserializer)
    : rho(T())
{
    this->unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
NoDynamics<T,Descriptor>* NoDynamics<T,Descriptor>::clone() const {
    return new NoDynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int NoDynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    Dynamics<T,Descriptor>::serialize(serializer);
    serializer.addValue(rho);
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    PLB_PRECONDITION( unserializer.getId() == this->getId() );
    Dynamics<T,Descriptor>::unserialize(unserializer);
    rho = unserializer.readValue<T>();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{ }

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar,
        Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat )
{ }

template<typename T, template<typename U> class Descriptor>
T NoDynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                            T jSqr, T thetaBar) const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::regularize(
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar ) const
{ }


template<typename T, template<typename U> class Descriptor>
T NoDynamics<T,Descriptor>::computeDensity(Cell<T,Descriptor> const& cell) const {
    return rho;
}

template<typename T, template<typename U> class Descriptor>
T NoDynamics<T,Descriptor>::computePressure(Cell<T,Descriptor> const& cell) const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::computeVelocity (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& u) const
{
    u.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
T NoDynamics<T,Descriptor>::computeTemperature(Cell<T,Descriptor> const& cell) const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::computePiNeq (
        Cell<T,Descriptor> const& cell,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    PiNeq.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::computeShearStress (
        Cell<T,Descriptor> const& cell,
        Array<T,SymmetricTensor<T,Descriptor>::n>& stress ) const
{
    stress.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::computeHeatFlux (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& q) const
{
    q.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::computeMoment (
        Cell<T,Descriptor> const& cell, plint momentId, T* moment) const
{ }

template<typename T, template<typename U> class Descriptor>
T NoDynamics<T,Descriptor>::getOmega() const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::setOmega(T omega_)
{ }

template<typename T, template<typename U> class Descriptor>
T NoDynamics<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const
{
    return (T)1 - Descriptor<T>::SkordosFactor();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j) const
{
    rhoBar = (T)1 - Descriptor<T>::SkordosFactor();
    j.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::computeRhoBarJPiNeq (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    rhoBar = (T)1 - Descriptor<T>::SkordosFactor();
    j.resetToZero();
    PiNeq.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
T NoDynamics<T,Descriptor>::computeEbar(Cell<T,Descriptor> const& cell) const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
plint NoDynamics<T,Descriptor>::numDecomposedVariables(plint order) const {
    return Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars;
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::decompose (
        Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order ) const
{
    rawData.resize(numDecomposedVariables(order));
    std::fill(rawData.begin(), rawData.end(), T());
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::recompose (
        Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order ) const
{
    PLB_PRECONDITION( (plint)rawData.size() == numDecomposedVariables(order) );
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::rescale (
        std::vector<T>& rawData, T xDxInv, T xDt, plint order ) const
{ }


template<typename T, template<typename U> class Descriptor>
bool NoDynamics<T,Descriptor>::hasMoments() const
{
    return false;
}


/******************************************************************************/

template<typename T, template<typename U> class Descriptor>
void constructIdChain(Dynamics<T,Descriptor> const& dynamics, std::vector<int>& chain)
{
    chain.push_back(dynamics.getId());
    if(dynamics.isComposite()) {
        constructIdChain(dynamic_cast<CompositeDynamics<T,Descriptor> const&>(dynamics).getBaseDynamics(), chain);
    }
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor> const& getBottomMostDynamics(Dynamics<T,Descriptor> const& dynamics)
{
    if(dynamics.isComposite()) {
        return getBottomMostDynamics (
                dynamic_cast<CompositeDynamics<T,Descriptor> const&>(dynamics).getBaseDynamics() );
    }
    else {
        return dynamics;
    }
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>* cloneAndReplaceBottomDynamics(Dynamics<T,Descriptor> const& dynamics,
                                                      Dynamics<T,Descriptor>* newBottom)
{
    if (!dynamics.isComposite()) {
        return newBottom;
    }
    Dynamics<T,Descriptor>* clonedDynamics = dynamics.clone();
    Dynamics<T,Descriptor>* component = clonedDynamics;
    Dynamics<T,Descriptor>* child = clonedDynamics;
    while (child->isComposite()) {
        component = child;
        child = &(dynamic_cast<CompositeDynamics<T,Descriptor>*>(child)->getBaseDynamics());
    }
    dynamic_cast<CompositeDynamics<T,Descriptor>*>(component)->replaceBaseDynamics(newBottom);
    return clonedDynamics;
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>* cloneAndInsertAtTopDynamics(Dynamics<T,Descriptor> const& dynamics,
                                                    CompositeDynamics<T,Descriptor>* newTop)
{
    if (dynamics.isComposeable()) {
        newTop->replaceBaseDynamics(dynamics.clone());
    }
    else {
        // Original top dynamics is composite but not composeable.
        if (dynamics.isComposite()) {
            CompositeDynamics<T,Descriptor> const& compositeDynamics =
                dynamic_cast<CompositeDynamics<T,Descriptor>const&>(dynamics);
            newTop->replaceBaseDynamics(compositeDynamics.getBaseDynamics().clone());
            // New top is composeable. Insert right behind original top dynamics.
            if (newTop->isComposeable()) {
                newTop = compositeDynamics.cloneWithNewBase(newTop);
            }
            // If new top is not composeable, then ignore original top dynamics.
        }
        // If original top dynamics is non composite and not composeable,
        //   then ignore original top dynamics.
    }
    return newTop;
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>* removeBoundaryComponents(Dynamics<T,Descriptor> const& dynamics)
{
    if (dynamics.isComposite()) {
        CompositeDynamics<T,Descriptor> const& compositeDynamics =
            dynamic_cast<CompositeDynamics<T,Descriptor>const&>(dynamics);
        Dynamics<T,Descriptor>* noBoundaryBase = removeBoundaryComponents(compositeDynamics.getBaseDynamics());
        if (dynamics.isBoundary()) {
            return noBoundaryBase;
        }
        else {
            return  dynamics.cloneWithNewBase(noBoundaryBase);
        }
    }
    else {
        /// At bottom-most level, dynamics is kept even if it is boundary (e.g. Bounce-Back).
        return dynamics.clone();
    }
}

template<typename T, template<typename U> class Descriptor>
void serialize(Dynamics<T,Descriptor> const& dynamics, std::vector<char>& data) {
    HierarchicSerializer serializer(data, dynamics.getId());
    dynamics.serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void serialize(std::vector<Dynamics<T,Descriptor>*> const& dynamics, std::vector<char>& data)
{
    for (pluint iDyn=0; iDyn<dynamics.size(); ++iDyn) {
        serialize(*dynamics[iDyn], data);
    }
}

template<typename T, template<typename U> class Descriptor>
pluint unserialize( Dynamics<T,Descriptor>& dynamics,
                    std::vector<char> const& data, pluint serializerPos )
{
    PLB_ASSERT( serializerPos < data.size() );
    HierarchicUnserializer unserializer(data, serializerPos);
    dynamics.unserialize(unserializer);
    return unserializer.getCurrentPos();
}

/// Unserialize all data into newly generated dynamics objects.
template<typename T, template<typename U> class Descriptor>
void generateAndUnserializeDynamics (
        std::vector<char> const& data,
        std::vector<Dynamics<T,Descriptor>*>& dynamics)
{
    pluint serializerPos = 0;
    while (serializerPos < data.size()) {
        HierarchicUnserializer unserializer(data, serializerPos);
        dynamics.push_back (
                meta::dynamicsRegistration<T,Descriptor>().generate(unserializer) );
        serializerPos = unserializer.getCurrentPos();
    }
}


}  // namespace plb

#endif  // DYNAMICS_HH
