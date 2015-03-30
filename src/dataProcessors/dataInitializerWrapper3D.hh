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

/** \file
 * Helper functions for domain initialization -- generic implementation.
 */
#ifndef DATA_INITIALIZER_WRAPPER_3D_HH
#define DATA_INITIALIZER_WRAPPER_3D_HH

#include "dataProcessors/dataInitializerWrapper3D.h"
#include "core/cell.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Initialization of the block-lattice: atomic-block * */
/* ******************************************************************* */

template<typename T, template<class U> class Descriptor>
void apply(BlockLattice3D<T,Descriptor>& lattice, Box3D domain, OneCellFunctional3D<T,Descriptor>* f) {
    applyProcessingFunctional(new GenericLatticeFunctional3D<T,Descriptor>(f), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void applyIndexed(BlockLattice3D<T,Descriptor>& lattice, Box3D domain, OneCellIndexedFunctional3D<T,Descriptor>* f) {
    applyProcessingFunctional(new GenericIndexedLatticeFunctional3D<T,Descriptor>(f), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLattice3D<T,Descriptor>& lattice, Box3D domain, Dynamics<T,Descriptor>* dynamics) {
    applyProcessingFunctional (
        new InstantiateDynamicsFunctional3D<T,Descriptor>(dynamics), domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLattice3D<T,Descriptor>& lattice, Box3D boundingBox,
                    DomainFunctional3D* domain, Dynamics<T,Descriptor>* dynamics) {
    applyProcessingFunctional (
        new InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>(dynamics, domain),
        boundingBox, lattice );
}

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLattice3D<T,Descriptor>& lattice,
                    DotList3D const& dotList, Dynamics<T,Descriptor>* dynamics)
{
    applyProcessingFunctional (
        new InstantiateDotDynamicsFunctional3D<T,Descriptor>(dynamics), dotList, lattice );
}

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLattice3D<T,Descriptor>& lattice,
                    plint iX, plint iY, plint iZ, Dynamics<T,Descriptor>* dynamics)
{
    DotList3D pos; pos.addDot(Dot3D(iX,iY,iZ));
    defineDynamics(lattice, pos, dynamics);
}

template<typename T, template<typename U> class Descriptor>
void defineDynamics( BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<bool>& boolMask,
                     Box3D domain, Dynamics<T,Descriptor>* dynamics, bool whichFlag )
{
    applyProcessingFunctional (
            new DynamicsFromMaskFunctional3D<T,Descriptor>(dynamics, whichFlag),
            domain, lattice, boolMask );
}

template<typename T, template<typename U> class Descriptor>
void defineDynamics( BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<bool>& boolMask,
                     Dynamics<T,Descriptor>* dynamics, bool whichFlag )
{
    defineDynamics(lattice, boolMask, lattice.getBoundingBox(), dynamics, whichFlag);
}


template<typename T, template<typename U> class Descriptor>
void defineDynamics( BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<int>& intMask,
                     Box3D domain, Dynamics<T,Descriptor>* dynamics, int whichFlag )
{
    applyProcessingFunctional (
            new DynamicsFromIntMaskFunctional3D<T,Descriptor>(dynamics, whichFlag),
            domain, lattice, intMask );
}

template<typename T, template<typename U> class Descriptor>
void defineDynamics( BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<int>& intMask,
                     Dynamics<T,Descriptor>* dynamics, int whichFlag )
{
    defineDynamics(lattice, intMask, lattice.getBoundingBox(), dynamics, whichFlag);
}

template<typename T, template<typename U> class Descriptor>
void recomposeFromFlowVariables ( BlockLattice3D<T,Descriptor>& lattice,
                                  ScalarField3D<T>& density, TensorField3D<T,3>& velocity,
                                  TensorField3D<T,6>& strainRate, Box3D domain )
{
    std::vector<AtomicBlock3D*> atomicBlocks(4);
    atomicBlocks[0] = &lattice;
    atomicBlocks[1] = &density;
    atomicBlocks[2] = &velocity;
    atomicBlocks[3] = &strainRate;
    applyProcessingFunctional( new RecomposeFromFlowVariablesFunctional3D<T,Descriptor>, domain,
                               atomicBlocks );
}

template<typename T, template<typename U> class Descriptor>
void recomposeFromFlowVariables ( BlockLattice3D<T,Descriptor>& lattice,
                                  ScalarField3D<T>& density, TensorField3D<T,3>& velocity,
                                  TensorField3D<T,6>& strainRate )
{
    recomposeFromFlowVariables( lattice, density, velocity, strainRate, lattice.getBoundingBox() );
}

template<typename T, template<class U> class Descriptor>
void setOmega(BlockLattice3D<T,Descriptor>& lattice, Box3D domain, T omega) {
    applyProcessingFunctional(new AssignOmegaFunctional3D<T,Descriptor>(omega), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void setOmega(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T> &omega, Box3D domain) {
    applyProcessingFunctional(new AssignScalarFieldOmegaFunctional3D<T,Descriptor>(), domain, lattice, omega);
}

template<typename T, template<class U> class Descriptor>
void setBoundaryVelocity(BlockLattice3D<T,Descriptor>& lattice, Box3D domain, Array<T,3> velocity) {
    applyProcessingFunctional(new SetConstBoundaryVelocityFunctional3D<T,Descriptor>(velocity), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void setBoundaryDensity(BlockLattice3D<T,Descriptor>& lattice, Box3D domain, T rho) {
    applyProcessingFunctional(new SetConstBoundaryDensityFunctional3D<T,Descriptor>(rho), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void initializeAtEquilibrium(BlockLattice3D<T,Descriptor>& lattice, Box3D domain, T rho, Array<T,3> velocity) {
    applyProcessingFunctional(new IniConstEquilibriumFunctional3D<T,Descriptor>(rho, velocity), domain, lattice);
}

template<typename T, template<class U> class Descriptor, class DomainFunctional>
void maskedInitializeAtEquilibrium(MultiBlockLattice3D<T,Descriptor>& lattice,
                                   Box3D boundingBox, DomainFunctional const& domain,
                                   T density, Array<T,3> velocity)
{
    applyProcessingFunctional(new IniConstEquilibriumOnDomainFunctional3D<T,Descriptor,DomainFunctional>(density, velocity, (T)1, domain), boundingBox, lattice);
}

template<typename T, template<class U> class Descriptor>
void stripeOffDensityOffset(BlockLattice3D<T,Descriptor>& lattice, Box3D domain, T deltaRho) {
    applyProcessingFunctional(new StripeOffDensityOffsetFunctional3D<T,Descriptor>(deltaRho), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void setCompositeDynamics( BlockLattice3D<T,Descriptor>& lattice, Box3D domain,
                           CompositeDynamics<T,Descriptor>* compositeDynamics )
{
    applyProcessingFunctional (
            new InstantiateCompositeDynamicsFunctional3D<T,Descriptor>(compositeDynamics),
            domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void setExternalScalar( BlockLattice3D<T,Descriptor>& lattice, Box3D domain,
                        int whichScalar, T externalScalar )
{
    applyProcessingFunctional (
            new SetExternalScalarFunctional3D<T,Descriptor>(whichScalar, externalScalar),
            domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void setExternalScalar( BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<int>& mask,
                        int flag, Box3D domain, int whichScalar, T externalScalar )
{
    applyProcessingFunctional (
            new MaskedSetExternalScalarFunctional3D<T,Descriptor>(flag, whichScalar, externalScalar),
            domain, lattice, mask );
}

template<typename T, template<class U> class Descriptor>
void setExternalVector( BlockLattice3D<T,Descriptor>& lattice, Box3D domain,
                        int vectorStartsAt, Array<T,Descriptor<T>::d> externalVector )
{
    applyProcessingFunctional (
            new SetExternalVectorFunctional3D<T,Descriptor>(vectorStartsAt, externalVector),
            domain, lattice );
}



/* *************** PART II ******************************************* */
/* *************** Initialization of the block-lattice: multi-block * */
/* ******************************************************************* */

template<typename T, template<class U> class Descriptor>
void apply(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain, OneCellFunctional3D<T,Descriptor>* f) {
    applyProcessingFunctional(new GenericLatticeFunctional3D<T,Descriptor>(f), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void applyIndexed(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain, OneCellIndexedFunctional3D<T,Descriptor>* f) {
    applyProcessingFunctional(new GenericIndexedLatticeFunctional3D<T,Descriptor>(f), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void defineDynamics(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain, Dynamics<T,Descriptor>* dynamics) {
    applyProcessingFunctional (
        new InstantiateDynamicsFunctional3D<T,Descriptor>(dynamics), domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void defineDynamics(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D boundingBox,
                    DomainFunctional3D* domain, Dynamics<T,Descriptor>* dynamics) {
    applyProcessingFunctional (
        new InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>(dynamics, domain),
        boundingBox, lattice );
}

template<typename T, template<class U> class Descriptor>
void defineDynamics(MultiBlockLattice3D<T,Descriptor>& lattice,
                    DotList3D const& dotList, Dynamics<T,Descriptor>* dynamics)
{
    applyProcessingFunctional (
        new InstantiateDotDynamicsFunctional3D<T,Descriptor>(dynamics), dotList, lattice );
}

template<typename T, template<class U> class Descriptor>
void defineDynamics(MultiBlockLattice3D<T,Descriptor>& lattice,
                    plint iX, plint iY, plint iZ, Dynamics<T,Descriptor>* dynamics)
{
    DotList3D pos; pos.addDot(Dot3D(iX,iY,iZ));
    defineDynamics(lattice, pos, dynamics);
}

template<typename T, template<typename U> class Descriptor>
void defineDynamics( MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<bool>& boolMask,
                     Box3D domain, Dynamics<T,Descriptor>* dynamics, bool whichFlag )
{
    applyProcessingFunctional (
            new DynamicsFromMaskFunctional3D<T,Descriptor>(dynamics, whichFlag),
            domain, lattice, boolMask );
}

template<typename T, template<typename U> class Descriptor>
void defineDynamics( MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<bool>& boolMask,
                     Dynamics<T,Descriptor>* dynamics, bool whichFlag )
{
    defineDynamics(lattice, boolMask, lattice.getBoundingBox(), dynamics, whichFlag);
}

template<typename T, template<typename U> class Descriptor>
void defineDynamics( MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<int>& intMask,
                     Box3D domain, Dynamics<T,Descriptor>* dynamics, int whichFlag )
{
    applyProcessingFunctional (
            new DynamicsFromIntMaskFunctional3D<T,Descriptor>(dynamics, whichFlag),
            domain, lattice, intMask );
}

template<typename T, template<typename U> class Descriptor>
void defineDynamics( MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<int>& intMask,
                     Dynamics<T,Descriptor>* dynamics, int whichFlag )
{
    defineDynamics(lattice, intMask, lattice.getBoundingBox(), dynamics, whichFlag);
}

template<typename T, template<typename U> class Descriptor>
void recomposeFromFlowVariables ( MultiBlockLattice3D<T,Descriptor>& lattice,
                                  MultiScalarField3D<T>& density, MultiTensorField3D<T,3>& velocity,
                                  MultiTensorField3D<T,6>& strainRate, Box3D domain )
{
    std::vector<MultiBlock3D*> multiBlocks(4);
    multiBlocks[0] = &lattice;
    multiBlocks[1] = &density;
    multiBlocks[2] = &velocity;
    multiBlocks[3] = &strainRate;
    applyProcessingFunctional( new RecomposeFromFlowVariablesFunctional3D<T,Descriptor>, domain,
                               multiBlocks );
}

template<typename T, template<typename U> class Descriptor>
void recomposeFromOrderZeroVariables ( MultiBlockLattice3D<T,Descriptor>& lattice,
                                  MultiScalarField3D<T>& density, MultiTensorField3D<T,3>& velocity,
                                  MultiTensorField3D<T,Descriptor<T>::q>& fNeq )
{
    recomposeFromOrderZeroVariables( lattice, density, velocity, fNeq, lattice.getBoundingBox() );
}


template<typename T, template<typename U> class Descriptor>
void recomposeFromOrderZeroVariables ( MultiBlockLattice3D<T,Descriptor>& lattice,
                                  MultiScalarField3D<T>& density, MultiTensorField3D<T,3>& velocity,
                                  MultiTensorField3D<T,Descriptor<T>::q>& fNeq, Box3D domain )
{
    std::vector<MultiBlock3D*> multiBlocks(4);
    multiBlocks[0] = &lattice;
    multiBlocks[1] = &density;
    multiBlocks[2] = &velocity;
    multiBlocks[3] = &fNeq;
    applyProcessingFunctional( new RecomposeFromOrderZeroVariablesFunctional3D<T,Descriptor>, domain,
                               multiBlocks );
}

template<typename T, template<typename U> class Descriptor>
void recomposeFromFlowVariables ( MultiBlockLattice3D<T,Descriptor>& lattice,
                                  MultiScalarField3D<T>& density, MultiTensorField3D<T,3>& velocity,
                                  MultiTensorField3D<T,6>& strainRate )
{
    recomposeFromFlowVariables( lattice, density, velocity, strainRate, lattice.getBoundingBox() );
}

template<typename T, template<class U> class Descriptor>
void setOmega(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain, T omega) {
    applyProcessingFunctional(new AssignOmegaFunctional3D<T,Descriptor>(omega), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void setOmega(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& omega, Box3D domain) {
    applyProcessingFunctional(new AssignScalarFieldOmegaFunctional3D<T,Descriptor>(), domain, lattice, omega);
}

template<typename T, template<class U> class Descriptor>
void setBoundaryVelocity(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain, Array<T,3> velocity) {
    applyProcessingFunctional(new SetConstBoundaryVelocityFunctional3D<T,Descriptor>(velocity), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void setBoundaryDensity(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain, T rho) {
    applyProcessingFunctional(new SetConstBoundaryDensityFunctional3D<T,Descriptor>(rho), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void initializeAtEquilibrium(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain, T rho, Array<T,3> velocity) {
    applyProcessingFunctional(new IniConstEquilibriumFunctional3D<T,Descriptor>(rho, velocity, (T)1), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void initializeAtThermalEquilibrium(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain, 
                                    T rho, Array<T,3> velocity, T temperature) {
    applyProcessingFunctional(new IniConstEquilibriumFunctional3D<T,Descriptor>(rho, velocity, temperature), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void stripeOffDensityOffset(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain, T deltaRho) {
    applyProcessingFunctional(new StripeOffDensityOffsetFunctional3D<T,Descriptor>(deltaRho), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void setCompositeDynamics( MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
                           CompositeDynamics<T,Descriptor>* compositeDynamics )
{
    applyProcessingFunctional (
            new InstantiateCompositeDynamicsFunctional3D<T,Descriptor>(compositeDynamics),
            domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void setExternalScalar( MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
                        int whichScalar, T externalScalar )
{
    applyProcessingFunctional (
            new SetExternalScalarFunctional3D<T,Descriptor>(whichScalar, externalScalar),
            domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void setExternalScalar( MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
                        int whichScalar, MultiScalarField3D<T> &scalar )
{
    applyProcessingFunctional (
            new SetExternalScalarFromScalarFieldFunctional3D<T,Descriptor>(whichScalar),
            domain, lattice, scalar );
}

template<typename T, template<class U> class Descriptor>
void setExternalScalar( MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<int>& mask,
                        int flag, Box3D domain, int whichScalar, T externalScalar )
{
    applyProcessingFunctional (
            new MaskedSetExternalScalarFunctional3D<T,Descriptor>(flag, whichScalar, externalScalar),
            domain, lattice, mask );
}

template<typename T, template<class U> class Descriptor, class Functional>
void setGenericExternalScalar( MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
                               int whichScalar, Functional const& functional )
{
    applyProcessingFunctional (
            new SetGenericExternalScalarFunctional3D<T,Descriptor,Functional>(whichScalar, functional),
            domain, lattice );
}

template<typename T, template<class U> class Descriptor, class Functional>
void setGenericExternalScalar( MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<int>& mask,
                               int flag, Box3D domain, int whichScalar, Functional const& functional )
{
    applyProcessingFunctional (
            new MaskedSetGenericExternalScalarFunctional3D<T,Descriptor,Functional>(flag, whichScalar, functional),
            domain, lattice, mask );
}

template<typename T, template<class U> class Descriptor>
void setExternalVector( MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
                        int vectorStartsAt, Array<T,Descriptor<T>::d> externalVector )
{
    applyProcessingFunctional (
            new SetExternalVectorFunctional3D<T,Descriptor>(vectorStartsAt, externalVector),
            domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void setExternalVector( MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<int>& mask,
                        int flag, Box3D domain, int vectorStartsAt, Array<T,Descriptor<T>::d> externalVector )
{
    applyProcessingFunctional (
            new MaskedSetExternalVectorFunctional3D<T,Descriptor>(flag, vectorStartsAt, externalVector),
            domain, lattice, mask );
}

template<typename T, template<class U> class Descriptor, int nDim>
void setExternalVector( MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
                        int vectorStartsAt, MultiTensorField3D<T,nDim> &tensor )
{
    applyProcessingFunctional (
            new SetExternalVectorFromTensorFieldFunctional3D<T,Descriptor,nDim>(vectorStartsAt),
            domain, lattice, tensor );
}

template<typename T, template<class U> class Descriptor, class Functional>
void setExternalVector( MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
                        int vectorStartsAt, Functional const& functional )
{
    applyProcessingFunctional (
        new SetGenericExternalVectorFunctional3D<T,Descriptor,Functional>(vectorStartsAt, functional),
                               domain, lattice );
}

template<typename T, template<class U> class Descriptor, class Functional>
void setExternalVector( MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<int>& mask,
                        int flag, Box3D domain, int vectorStartsAt, Functional const& functional )
{
    applyProcessingFunctional (
        new MaskedSetGenericExternalVectorFunctional3D<T,Descriptor,Functional>(flag, vectorStartsAt, functional),
                               domain, lattice, mask );
}

template<typename T, template<class U> class Descriptor>
void interpolatePopulations (
        MultiBlockLattice3D<T,Descriptor>& modifiedLattice,
        MultiBlockLattice3D<T,Descriptor>& constLattice,
        plint minIter, plint maxIter, Box3D const& domain )
{
    applyProcessingFunctional (
            new InterpolatePopulationsFunctional3D<T,Descriptor>(minIter,maxIter),
            domain, modifiedLattice, constLattice );
}

template<typename T, template<class U> class Descriptor>
void interpolatePopulations (
        MultiBlockLattice3D<T,Descriptor>& modifiedLattice,
        MultiBlockLattice3D<T,Descriptor>& constLattice,
        plint minIter, plint maxIter )
{
    interpolatePopulations( modifiedLattice, constLattice, minIter, maxIter,
                            modifiedLattice.getBoundingBox() );
}

/* *************** PART III ******************************************* */
/* *************** Initialization of scalar- and tensor-fields: ****** */
/* *************** Atomic-Block  ************************************* */
/* ******************************************************************* */

template<typename T>
void setToConstant(ScalarField3D<T>& field, Box3D domain, T value) {
    applyProcessingFunctional(new IniConstScalarFunctional3D<T>(value), domain, field);
}

template<typename T>
void setToConstant( ScalarField3D<T>& field, ScalarField3D<int>& mask,
                    int flag, Box3D domain, T value )
{
    applyProcessingFunctional (
            new MaskedIniConstScalarFunctional3D<T>(flag, value), domain, field, mask);
}

template<typename T, int nDim>
void setToConstant( TensorField3D<T,nDim>& field, Box3D domain,
                    Array<T,nDim> const& value )
{
    applyProcessingFunctional (
            new IniConstTensorFunctional3D<T,nDim>(value), domain, field );
}

template<typename T, int nDim>
void setToConstant( TensorField3D<T,nDim>& field, ScalarField3D<int>& mask, int flag,
                    Box3D domain, Array<T,nDim> const& value )
{
    applyProcessingFunctional (
            new MaskedIniConstTensorFunctional3D<T,nDim>(flag, value), domain, mask, field );
}

template<typename T>
void setToCoordinate(ScalarField3D<T>& field, Box3D domain, plint index) {
    applyProcessingFunctional(new SetToCoordinateFunctional3D<T>(index), domain, field);
}

template<typename T>
void setToCoordinates(TensorField3D<T,3>& field, Box3D domain) {
    applyProcessingFunctional(new SetToCoordinatesFunctional3D<T>, domain, field);
}

template<typename T, int nDim>
void assignComponent(TensorField3D<T,nDim>& tensorField, int whichComponent,
                     ScalarField3D<T>& scalarField, Box3D domain)
{
    applyProcessingFunctional(new SetTensorComponentFunctional3D<T,nDim>(whichComponent),
                              domain, scalarField, tensorField);
}

/* *************** PART IV ******************************************* */
/* *************** Initialization of scalar- and tensor-fields: ****** */
/* *************** Multi-Block  ************************************** */
/* ******************************************************************* */

template<typename T>
void setToConstant(MultiScalarField3D<T>& field, Box3D domain, T value) {
    applyProcessingFunctional(new IniConstScalarFunctional3D<T>(value), domain, field);
}

template<typename T>
void setToConstant( MultiScalarField3D<T>& field, MultiScalarField3D<int>& mask,
                    int flag, Box3D domain, T value )
{
    applyProcessingFunctional (
            new MaskedIniConstScalarFunctional3D<T>(flag, value), domain, field, mask);
}

template<typename T, int nDim>
void setToConstant( MultiTensorField3D<T,nDim>& field, Box3D domain,
                    Array<T,nDim> const& value )
{
    applyProcessingFunctional (
            new IniConstTensorFunctional3D<T,nDim>(value), domain, field );
}

template<typename T, int nDim>
void setToConstant( MultiTensorField3D<T,nDim>& field, MultiScalarField3D<int>& mask, int flag,
                    Box3D domain, Array<T,nDim> const& value )
{
    applyProcessingFunctional (
            new MaskedIniConstTensorFunctional3D<T,nDim>(flag, value), domain, mask, field );
}

template<typename T>
void setToCoordinate(MultiScalarField3D<T>& field, Box3D domain, plint index) {
    applyProcessingFunctional(new SetToCoordinateFunctional3D<T>(index), domain, field);
}

template<typename T>
void setToCoordinates(MultiTensorField3D<T,3>& field, Box3D domain) {
    applyProcessingFunctional(new SetToCoordinatesFunctional3D<T>, domain, field);
}

template<typename T, int nDim>
void assignComponent(MultiTensorField3D<T,nDim>& tensorField, int whichComponent,
                     MultiScalarField3D<T>& scalarField, Box3D domain)
{
    applyProcessingFunctional(new SetTensorComponentFunctional3D<T,nDim>(whichComponent),
                              domain, scalarField, tensorField);
}

}  // namespace plb

#endif  // DATA_INITIALIZER_WRAPPER_3D_HH
