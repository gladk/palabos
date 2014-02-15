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
 * Helper functions for data analysis -- header file.
 */
#ifndef DATA_ANALYSIS_WRAPPER_3D_HH
#define DATA_ANALYSIS_WRAPPER_3D_HH

#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "dataProcessors/dataAnalysisFunctional3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "multiBlock/multiBlockGenerator3D.h"

namespace plb {

/* ******************************************************************* */
/* *************** PART I. Atomic-block wrappers: Block-Lattice ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template<typename T, template<typename U> class Descriptor> 
T computeAverageDensity(BlockLattice3D<T,Descriptor>& lattice, Box3D domain) {
    BoxSumRhoBarFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return Descriptor<T>::fullRho( functional.getSumRhoBar() / (T) domain.nCells() );
}

template<typename T, template<typename U> class Descriptor> 
T computeAverageDensity(BlockLattice3D<T,Descriptor>& lattice) {
    return computeAverageDensity(lattice, lattice.getBoundingBox());
}


template<typename T, template<typename U> class Descriptor> 
T computeAverageRhoBar(BlockLattice3D<T,Descriptor>& lattice, Box3D domain) {
    BoxSumRhoBarFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumRhoBar() / (T) domain.nCells();
}

template<typename T, template<typename U> class Descriptor> 
T computeAverageRhoBar(BlockLattice3D<T,Descriptor>& lattice) {
    return computeAverageRhoBar(lattice, lattice.getBoundingBox());
}


template<typename T, template<typename U> class Descriptor> 
T computeAverageEnergy(BlockLattice3D<T,Descriptor>& lattice, Box3D domain) 
{
    BoxSumEnergyFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumEnergy() / (T) domain.nCells();;
}

template<typename T, template<typename U> class Descriptor> 
T computeAverageEnergy(BlockLattice3D<T,Descriptor>& lattice) {
    return computeAverageEnergy(lattice, lattice.getBoundingBox());
}

template<typename T, template<typename U> class Descriptor, class BoolMask> 
plint count(BlockLattice3D<T,Descriptor>& lattice, Box3D domain, BoolMask boolMask)
{
    CountLatticeElementsFunctional3D<T,Descriptor,BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getCount();
}

template<typename T, template<typename U> class Descriptor, class BoolMask> 
plint count(BlockLattice3D<T,Descriptor>& lattice, BoolMask boolMask)
{
    return count(lattice, lattice.getBoundingBox(), boolMask);
}

/* *************** Density ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeDensity(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& density)
{
    applyProcessingFunctional (
            new BoxDensityFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, density );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeDensity(BlockLattice3D<T,Descriptor>& lattice)
{
    ScalarField3D<T>* density
        = new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeDensity(lattice, *density);
    return std::auto_ptr<ScalarField3D<T> >(density);
}


/* *************** RhoBar ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeRhoBar(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& rhoBar)
{
    applyProcessingFunctional (
            new BoxRhoBarFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, rhoBar );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeRhoBar(BlockLattice3D<T,Descriptor>& lattice)
{
    ScalarField3D<T>* rhoBar
        = new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeRhoBar(lattice, *rhoBar);
    return std::auto_ptr<ScalarField3D<T> >(rhoBar);
}


/* *************** Kinetic Energy ************************************ */

template<typename T, template<typename U> class Descriptor>
void computeKineticEnergy(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& energy)
{
    applyProcessingFunctional (
            new BoxKineticEnergyFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, energy );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeKineticEnergy(BlockLattice3D<T,Descriptor>& lattice)
{
    ScalarField3D<T>* energy
        = new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeKineticEnergy(lattice, *energy);
    return std::auto_ptr<ScalarField3D<T> >(energy);
}


/* *************** Velocity Norm ************************************* */

template<typename T, template<typename U> class Descriptor>
void computeVelocityNorm(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& velocityNorm)
{
    applyProcessingFunctional (
            new BoxVelocityNormFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, velocityNorm );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeVelocityNorm(BlockLattice3D<T,Descriptor>& lattice)
{
    ScalarField3D<T>* velocityNorm
        = new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeVelocityNorm(lattice, *velocityNorm);
    return std::auto_ptr<ScalarField3D<T> >(velocityNorm);
}


/* *************** Velocity Component ******************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocityComponent(BlockLattice3D<T,Descriptor>& lattice,
                              ScalarField3D<T>& velocityComponent,
                              plint iComponent)
{
    applyProcessingFunctional (
            new BoxVelocityComponentFunctional3D<T,Descriptor>(iComponent),
            lattice.getBoundingBox(), lattice, velocityComponent );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeVelocityComponent (
        BlockLattice3D<T,Descriptor>& lattice, plint iComponent )
{
    ScalarField3D<T>* velocityComponent
        = new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeVelocityComponent(lattice, *velocityComponent, iComponent);
    return std::auto_ptr<ScalarField3D<T> >(velocityComponent);
}


/* *************** Velocity ****************************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocity(BlockLattice3D<T,Descriptor>& lattice,
                     TensorField3D<T,Descriptor<T>::d>& velocity)
{
    applyProcessingFunctional (
            new BoxVelocityFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, velocity );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField3D<T,Descriptor<T>::d> > computeVelocity(BlockLattice3D<T,Descriptor>& lattice)
{
    TensorField3D<T,Descriptor<T>::d>* velocity
        = new TensorField3D<T,Descriptor<T>::d>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeVelocity(lattice, *velocity);
    return std::auto_ptr<TensorField3D<T,Descriptor<T>::d> >(velocity);
}



/* *************** Pi Neq ********************************* */

template<typename T, template<typename U> class Descriptor>
void computePiNeq(BlockLattice3D<T,Descriptor>& lattice,
                  TensorField3D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq)
{
    applyProcessingFunctional (
            new BoxPiNeqFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, PiNeq );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField3D<T,SymmetricTensor<T,Descriptor>::n> > computePiNeq(BlockLattice3D<T,Descriptor>& lattice)
{
    TensorField3D<T,SymmetricTensor<T,Descriptor>::n>* PiNeq
        = new TensorField3D<T,SymmetricTensor<T,Descriptor>::n>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computePiNeq(lattice, *PiNeq);
    return std::auto_ptr<TensorField3D<T,SymmetricTensor<T,Descriptor>::n> >(PiNeq);
}

/* *************** Shear Stress ********************************* */

template<typename T, template<typename U> class Descriptor>
void computeShearStress(BlockLattice3D<T,Descriptor>& lattice,
                  TensorField3D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq)
{
    applyProcessingFunctional (
        new BoxShearStressFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, PiNeq );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField3D<T,SymmetricTensor<T,Descriptor>::n> > computeShearStress(BlockLattice3D<T,Descriptor>& lattice)
{
    TensorField3D<T,SymmetricTensor<T,Descriptor>::n>* PiNeq
        = new TensorField3D<T,SymmetricTensor<T,Descriptor>::n>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeShearStress(lattice, *PiNeq);
    return std::auto_ptr<TensorField3D<T,SymmetricTensor<T,Descriptor>::n> >(PiNeq);
}


/* *************** Strain Rate from Stress *************************** */

template<typename T, template<typename U> class Descriptor>
void computeStrainRateFromStress(BlockLattice3D<T,Descriptor>& lattice,
                             TensorField3D<T,SymmetricTensor<T,Descriptor>::n>& S)
{
    applyProcessingFunctional (
            new BoxStrainRateFromStressFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, S );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField3D<T,SymmetricTensor<T,Descriptor>::n> > computeStrainRateFromStress(BlockLattice3D<T,Descriptor>& lattice)
{
    TensorField3D<T,SymmetricTensor<T,Descriptor>::n>* S
        = new TensorField3D<T,SymmetricTensor<T,Descriptor>::n>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeStrainRateFromStress(lattice, *S);
    return std::auto_ptr<TensorField3D<T,SymmetricTensor<T,Descriptor>::n> >(S);
}



/* *************** Population *************************************** */

template<typename T, template<typename U> class Descriptor>
void computePopulation(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& population, plint iPop)
{
    applyProcessingFunctional (
            new BoxPopulationFunctional3D<T,Descriptor>(iPop), lattice.getBoundingBox(), lattice, population );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computePopulation(BlockLattice3D<T,Descriptor>& lattice, plint iPop)
{
    ScalarField3D<T>* population = new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computePopulation(lattice, *population, iPop);
    return std::auto_ptr<ScalarField3D<T> >(population);
}

template<typename T, template<typename U> class Descriptor>
void computeAllPopulations(BlockLattice3D<T,Descriptor>& lattice, 
                           TensorField3D<T,Descriptor<T>::q>& populations,
                           Box3D domain)
{
    applyProcessingFunctional (
            new BoxAllPopulationsFunctional3D<T,Descriptor>(), domain, lattice, populations );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeAllPopulations(BlockLattice3D<T,Descriptor>& lattice)
{
    TensorField3D<T,Descriptor<T>::q>* populations = 
            new TensorField3D<T,Descriptor<T>::q>(lattice);
    computeAllPopulations(lattice, *populations);
    return std::auto_ptr<TensorField3D<T,Descriptor<T>::q> >(populations);
}

template<typename T, template<typename U> class Descriptor>
void copyPopulations(BlockLattice3D<T,Descriptor>& latticeFrom, BlockLattice3D<T,Descriptor>& latticeTo,
                     Box3D domain)
{
    applyProcessingFunctional (
            new CopyPopulationsFunctional3D<T,Descriptor>(), domain, latticeFrom, latticeTo );
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
void copyConvertPopulations(MultiBlockLattice3D<T1,Descriptor1>& latticeFrom, MultiBlockLattice3D<T2,Descriptor2>& latticeTo,
                            Box3D domain)
{
    applyProcessingFunctional( new CopyConvertPopulationsFunctional3D<T1,Descriptor1,T2,Descriptor2>(), domain, latticeFrom, latticeTo );
}

template<typename T, template<typename U> class Descriptor>
void copyAll(BlockLattice3D<T,Descriptor>& latticeFrom,
             BlockLattice3D<T,Descriptor>& latticeTo, Box3D domain)
{
    applyProcessingFunctional (
            new LatticeCopyAllFunctional3D<T,Descriptor>(), domain, latticeFrom, latticeTo );
}

template<typename T, template<typename U> class Descriptor>
void copyRegenerate(BlockLattice3D<T,Descriptor>& latticeFrom,
                    BlockLattice3D<T,Descriptor>& latticeTo, Box3D domain)
{
    applyProcessingFunctional (
            new LatticeRegenerateFunctional3D<T,Descriptor>(), domain, latticeFrom, latticeTo );
}


/* *************** Omega ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeOmega(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& omega, Box3D domain)
{
    applyProcessingFunctional (
        new BoxOmegaFunctional3D<T,Descriptor>, domain, lattice, omega );
}

template<typename T, template<typename U> class Descriptor>
void computeOmega(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& omega)
{
    computeOmega(lattice, omega, lattice.getBoundingBox());
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeOmega(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > omega =
    generateMultiScalarField<T>(lattice, domain);
    computeOmega(lattice, *omega, domain);
    return omega;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeOmega(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computeOmega(lattice, lattice.getBoundingBox());
}


/* *************** KinematicViscosity ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeKinematicViscosity(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& omega, Box3D domain)
{
    applyProcessingFunctional (
        new BoxKinematicViscosityFunctional3D<T,Descriptor>, domain, lattice, omega );
}

template<typename T, template<typename U> class Descriptor>
void computeKinematicViscosity(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& omega)
{
    computeKinematicViscosity(lattice, omega, lattice.getBoundingBox());
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeKinematicViscosity(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > omega =
    generateMultiScalarField<T>(lattice, domain);
    computeKinematicViscosity(lattice, *omega, domain);
    return omega;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeKinematicViscosity(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computeKinematicViscosity(lattice, lattice.getBoundingBox());
}

/* *************** ExternalForce ****************************************** */

template<typename T, template<typename U> class Descriptor>
void computeExternalForce(MultiBlockLattice3D<T,Descriptor>& lattice,
                     MultiTensorField3D<T,Descriptor<T>::d>& force, Box3D domain)
{
    applyProcessingFunctional (
            new BoxExternalForceFunctional3D<T,Descriptor>, domain, lattice, force );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::d> > computeExternalForce(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::d> > force
        = generateMultiTensorField<T,Descriptor<T>::d>(lattice, domain);
    computeExternalForce(lattice, *force, domain);
    return force;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::d> >
    computeExternalForce(MultiBlockLattice3D<T,Descriptor>& lattice)
{
    return computeExternalForce(lattice, lattice.getBoundingBox());
}



/* *************** ExternalScalar ****************************************** */

template<typename T, template<typename U> class Descriptor>
void computeExternalScalar(MultiBlockLattice3D<T,Descriptor>& lattice,
                     MultiScalarField3D<T>& scalar, int whichScalar, Box3D domain)
{
    applyProcessingFunctional (
            new BoxExternalScalarFunctional3D<T,Descriptor>(whichScalar), domain, lattice, scalar );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeExternalScalar(MultiBlockLattice3D<T,Descriptor>& lattice, int whichScalar, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > scalar
        = generateMultiScalarField<T>(lattice, domain);
    computeExternalScalar(lattice, *scalar, whichScalar, domain);
    return scalar;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> >
    computeExternalScalar(MultiBlockLattice3D<T,Descriptor>& lattice, int whichScalar)
{
    return computeExternalScalar(lattice, whichScalar, lattice.getBoundingBox());
}


/* ******************************************************************* */
/* *************** PART II. Atomic-block wrappers: Scalar-Field ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template<typename T>
T computeSum(ScalarField3D<T>& scalarField, Box3D domain) {
    BoxScalarSumFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar();
}

template<typename T>
T computeSum(ScalarField3D<T>& scalarField) {
    return computeSum(scalarField, scalarField.getBoundingBox());
}

template<typename T>
T computeBoundedSum(ScalarField3D<T>& scalarField, Box3D domain) {
    BoundedBoxScalarSumFunctional3D<T> functional;
    plint envelopeWidth=1;
    applyProcessingFunctional(functional, domain, scalarField, envelopeWidth);
    return functional.getSumScalar();
}

template<typename T>
T computeBoundedSum(ScalarField3D<T>& scalarField) {
    return computeBoundedSum(scalarField, scalarField.getBoundingBox());
}

template<typename T>
T computeAverage(ScalarField3D<T>& scalarField, Box3D domain) {
    BoxScalarSumFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar() / (T) domain.nCells();
}

template<typename T>
T computeAverage(ScalarField3D<T>& scalarField) {
    return computeAverage(scalarField, scalarField.getBoundingBox());
}


template<typename T>
T computeAverage(ScalarField3D<T>& scalarField, ScalarField3D<int>& mask, int flag, Box3D domain)
{
    MaskedBoxScalarAverageFunctional3D<T> functional(flag);
    applyProcessingFunctional(functional, domain, scalarField, mask);
    return functional.getAverageScalar();
}

template<typename T>
T computeAverage(ScalarField3D<T>& scalarField, ScalarField3D<int>& mask, int flag) {
    return computeAverage(scalarField, mask, flag, scalarField.getBoundingBox());
}


template<typename T>
T computeMin(ScalarField3D<T>& scalarField, Box3D domain) {
    BoxScalarMinFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getMinScalar();
}

template<typename T>
T computeMin(ScalarField3D<T>& scalarField) {
    return computeMin(scalarField, scalarField.getBoundingBox());
}


template<typename T>
T computeMax(ScalarField3D<T>& scalarField, Box3D domain) {
    BoxScalarMaxFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getMaxScalar();
}

template<typename T>
T computeMax(ScalarField3D<T>& scalarField) {
    return computeMax(scalarField, scalarField.getBoundingBox());
}


template<typename T>
T computeBoundedAverage(ScalarField3D<T>& scalarField, Box3D domain) {
    BoundedBoxScalarSumFunctional3D<T> functional;
    plint envelopeWidth=1;
    applyProcessingFunctional(functional, domain, scalarField, envelopeWidth);
    return functional.getSumScalar() /
             (T) ( (domain.getNx()-1)*(domain.getNy()-1)*(domain.getNz()-1) );
}

template<typename T>
T computeBoundedAverage(ScalarField3D<T>& scalarField) {
    return computeBoundedAverage(scalarField, scalarField.getBoundingBox());
}


template<typename T, class BoolMask> 
plint count(ScalarField3D<T>& field, Box3D domain, BoolMask boolMask)
{
    CountScalarElementsFunctional3D<T,BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

template<typename T, class BoolMask> 
plint count(ScalarField3D<T>& field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

/* *************** ScalarField - Scalar operations *************** */

template<typename T>
void add(ScalarField3D<T>& field, T scalar, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new A_plus_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > add(ScalarField3D<T>& field, T scalar)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    add(field, scalar, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void add(T scalar, ScalarField3D<T>& field, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new A_plus_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > add(T scalar, ScalarField3D<T>& field)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    add(scalar, field, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void subtract(ScalarField3D<T>& field, T scalar, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new A_minus_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > subtract(ScalarField3D<T>& field, T scalar)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    subtract(field, scalar, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void subtract(T scalar, ScalarField3D<T>& field, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new Alpha_minus_A_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > subtract(T scalar, ScalarField3D<T>& field)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    subtract(scalar, field, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void multiply(ScalarField3D<T>& field, T scalar, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new A_times_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > multiply(ScalarField3D<T>& field, T scalar)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    multiply(field, scalar, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void multiply(T scalar, ScalarField3D<T>& field, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new A_times_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > multiply(T scalar, ScalarField3D<T>& field)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    multiply(scalar, field, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void divide(ScalarField3D<T>& field, T scalar, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new A_dividedBy_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > divide(ScalarField3D<T>& field, T scalar)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    divide(field, scalar, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void divide(T scalar, ScalarField3D<T>& field, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new Alpha_dividedBy_A_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > divide(T scalar, ScalarField3D<T>& field)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    divide(scalar, field, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}

/* *************** ScalarField - Scalar inplace operations *************** */

template<typename T>
void addInPlace(ScalarField3D<T>& field, T scalar) {
    applyProcessingFunctional (
            new A_plus_alpha_inplace_functional3D<T>(scalar), field.getBoundingBox(), field);
}

template<typename T>
void subtractInPlace(ScalarField3D<T>& field, T scalar) {
    applyProcessingFunctional (
            new A_minus_alpha_inplace_functional3D<T>(scalar), field.getBoundingBox(), field);
}

template<typename T>
void multiplyInPlace(ScalarField3D<T>& field, T scalar) {
    applyProcessingFunctional (
            new A_times_alpha_inplace_functional3D<T>(scalar), field.getBoundingBox(), field);
}

template<typename T>
void divideInPlace(ScalarField3D<T>& field, T scalar) {
    applyProcessingFunctional (
            new A_dividedBy_alpha_inplace_functional3D<T>(scalar), field.getBoundingBox(), field);
}


/* *************** ScalarField - ScalarField operations *************** */

template<typename T1, typename T2>
void copy(ScalarField3D<T1>& field, ScalarField3D<T2>& convertedField)
{
    applyProcessingFunctional (
            new CopyConvertScalarFunctional3D<T1,T2>, field.getBoundingBox(), field, convertedField );
}

template<typename T1, typename T2>
std::auto_ptr<ScalarField3D<T2> > copyConvert(ScalarField3D<T1>& field)
{
    ScalarField3D<T2>* convertedField = new ScalarField3D<T2>(field.getNx(), field.getNy(), field.getNz());
    plb::copy(field, *convertedField);
    return std::auto_ptr<ScalarField3D<T2> >(convertedField);
}


template<typename T>
void add(ScalarField3D<T>& A, ScalarField3D<T>& B, ScalarField3D<T>& result)
{
    std::vector<ScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_plus_B_functional3D<T>, A.getBoundingBox(), fields );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > add(ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    add(A, B, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void subtract(ScalarField3D<T>& A, ScalarField3D<T>& B, ScalarField3D<T>& result)
{
    std::vector<ScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_minus_B_functional3D<T>, A.getBoundingBox(), fields );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > subtract(ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    subtract(A, B, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void multiply(ScalarField3D<T>& A, ScalarField3D<T>& B, ScalarField3D<T>& result)
{
    std::vector<ScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_times_B_functional3D<T>, A.getBoundingBox(), fields );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > multiply(ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    multiply(A, B, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void divide(ScalarField3D<T>& A, ScalarField3D<T>& B, ScalarField3D<T>& result)
{
    std::vector<ScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_dividedBy_B_functional3D<T>, A.getBoundingBox(), fields );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > divide(ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    divide(A, B, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}

/* *************** ScalarField operations *************** */

template<typename T>
void computeSqrt(ScalarField3D<T>& A, ScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeScalarSqrtFunctional3D<T>, domain, A, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeSqrt(ScalarField3D<T>& A, Box3D domain)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    computeSqrt(A, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeSqrt(ScalarField3D<T>& A)
{
    return computeSqrt(A, A.getBoundingBox());
}



template<typename T>
void computeAbsoluteValue(ScalarField3D<T>& A, ScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
    new ComputeAbsoluteValueFunctional3D<T>, domain, A, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeAbsoluteValue(ScalarField3D<T>& A, Box3D domain)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    computeAbsoluteValue(A, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeAbsoluteValue(ScalarField3D<T>& A)
{
    return computeAbsoluteValue(A, A.getBoundingBox());
}


/* *************** ScalarField - ScalarField inplace operations *************** */

template<typename T>
void addInPlace(ScalarField3D<T>& A, ScalarField3D<T>& B) {
    applyProcessingFunctional (
            new A_plus_B_inplace_functional3D<T>, A.getBoundingBox(), A, B );
}

template<typename T>
void subtractInPlace(ScalarField3D<T>& A, ScalarField3D<T>& B) {
    applyProcessingFunctional (
            new A_minus_B_inplace_functional3D<T>, A.getBoundingBox(), A, B );
}

template<typename T>
void multiplyInPlace(ScalarField3D<T>& A, ScalarField3D<T>& B) {
    applyProcessingFunctional (
            new A_times_B_inplace_functional3D<T>, A.getBoundingBox(), A, B );
}

template<typename T>
void divideInPlace(ScalarField3D<T>& A, ScalarField3D<T>& B) {
    applyProcessingFunctional (
            new A_dividedBy_B_inplace_functional3D<T>, A.getBoundingBox(), A, B );
}



/* ******************************************************************* */
/* *************** PART III. Atomic-block wrappers: Tensor-Field ***** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template<typename T, plint nDim, class BoolMask> 
plint count(TensorField3D<T,nDim>& field, Box3D domain, BoolMask boolMask)
{
    CountTensorElementsFunctional3D<T,nDim,BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

template<typename T, plint nDim, class BoolMask> 
plint count(TensorField3D<T,nDim>& field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

/* *************** Component (scalar-field) out of a tensor-field ****** */

template<typename T, int nDim>
void extractComponent(TensorField3D<T,nDim>& tensorField, ScalarField3D<T>& component, int iComponent)
{
    applyProcessingFunctional (
            new ExtractTensorComponentFunctional3D<T,nDim>(iComponent), tensorField.getBoundingBox(), component, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<ScalarField3D<T> > extractComponent(TensorField3D<T,nDim>& tensorField, int iComponent)
{
    ScalarField3D<T>* component = new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    extractComponent(tensorField, *component, iComponent);
    return std::auto_ptr<ScalarField3D<T> >(component);
}


/* *************** Vector-norm of each cell in the field *************** */

template<typename T, int nDim>
void computeNorm(TensorField3D<T,nDim>& tensorField, ScalarField3D<T>& component)
{
    applyProcessingFunctional (
            new ComputeNormFunctional3D<T,nDim>, tensorField.getBoundingBox(), component, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<ScalarField3D<T> > computeNorm(TensorField3D<T,nDim>& tensorField)
{
    ScalarField3D<T>* component = new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeNorm(tensorField, *component);
    return std::auto_ptr<ScalarField3D<T> >(component);
}

/* *************** Sqrt operation on each component of each cell *************** */

template<typename T, int nDim>
void computeSqrt(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& result, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeTensorSqrtFunctional3D<T, nDim>, domain, A, result );
}

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > computeSqrt(TensorField3D<T,nDim>& A, Box3D domain)
{
    TensorField3D<T,nDim>* result = new TensorField3D<T,nDim>(A.getNx(), A.getNy(), A.getNz());
    computeSqrt(A, *result, domain);
    return result;
}

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > computeSqrt(TensorField3D<T,nDim>& A)
{
    return computeSqrt(A, A.getBoundingBox());
}


/* *************** Squared vector-norm of each cell in the field ******** */

template<typename T, int nDim>
void computeNormSqr(TensorField3D<T,nDim>& tensorField, ScalarField3D<T>& component)
{
    applyProcessingFunctional (
            new ComputeNormSqrFunctional3D<T,nDim>, tensorField.getBoundingBox(), component, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<ScalarField3D<T> > computeNormSqr(TensorField3D<T,nDim>& tensorField)
{
    ScalarField3D<T>* component = new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeNormSqr(tensorField, *component);
    return std::auto_ptr<ScalarField3D<T> >(component);
}


/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template<typename T>
void computeSymmetricTensorNorm(TensorField3D<T,6>& tensorField, ScalarField3D<T>& norm)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorNormFunctional3D<T>, tensorField.getBoundingBox(), norm, tensorField );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeSymmetricTensorNorm(TensorField3D<T,6>& tensorField)
{
    ScalarField3D<T>* norm = new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeSymmetricTensorNorm(tensorField, *norm);
    return std::auto_ptr<ScalarField3D<T> >(norm);
}


/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template<typename T>
void computeSymmetricTensorNormSqr(TensorField3D<T,6>& tensorField, ScalarField3D<T>& normSqr)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorNormSqrFunctional3D<T>, tensorField.getBoundingBox(), normSqr, tensorField );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeSymmetricTensorNormSqr(TensorField3D<T,6>& tensorField)
{
    ScalarField3D<T>* normSqr = new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeSymmetricTensorNormSqr(tensorField, *normSqr);
    return std::auto_ptr<ScalarField3D<T> >(normSqr);
}


/* *************** Trace of each symmetric tensor of a field ************ */

template<typename T>
void computeSymmetricTensorTrace(TensorField3D<T,6>& tensorField, ScalarField3D<T>& trace)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorTraceFunctional3D<T>, tensorField.getBoundingBox(), trace, tensorField );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeSymmetricTensorTrace(TensorField3D<T,6>& tensorField)
{
    ScalarField3D<T>* trace = new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeSymmetricTensorTrace(tensorField, *trace);
    return std::auto_ptr<ScalarField3D<T> >(trace);
}




/* *************** Vorticity from Velocity field *********************** */

template<typename T>
void computeVorticity(TensorField3D<T,3>& velocity, TensorField3D<T,3>& vorticity)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxVorticityFunctional3D<T,3>, velocity.getBoundingBox(), velocity, vorticity, envelopeWidth );
}

template<typename T>
std::auto_ptr<TensorField3D<T,3> > computeVorticity(TensorField3D<T,3>& velocity)
{
    TensorField3D<T,3>* vorticity = new TensorField3D<T,3>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeVorticity(velocity, *vorticity);
    return std::auto_ptr<TensorField3D<T,3> >(vorticity);
}


/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkVorticity(TensorField3D<T,3>& velocity, TensorField3D<T,3>& vorticity)
{
    applyProcessingFunctional (
            new BoxBulkVorticityFunctional3D<T,3>, velocity.getBoundingBox(), velocity, vorticity );
}

template<typename T>
std::auto_ptr<TensorField3D<T,3> > computeBulkVorticity(TensorField3D<T,3>& velocity)
{
    TensorField3D<T,3>* vorticity = new TensorField3D<T,3>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeBulkVorticity(velocity, *vorticity);
    return std::auto_ptr<TensorField3D<T,3> >(vorticity);
}


/* *************** Divergence, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkDivergence(TensorField3D<T,3>& velocity, ScalarField3D<T>& divergence)
{
    applyProcessingFunctional (
            new BoxBulkDivergenceFunctional3D<T,3>, velocity.getBoundingBox(), divergence, velocity );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeBulkDivergence(TensorField3D<T,3>& velocity)
{
    ScalarField3D<T>* divergence = new ScalarField3D<T>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeBulkDivergence(velocity, *divergence);
    return std::auto_ptr<ScalarField3D<T> >(divergence);
}


/* *************** Strain rate from Velocity field ********************* */

template<typename T>
void computeStrainRate(TensorField3D<T,3>& velocity, TensorField3D<T,6>& S)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxStrainRateFunctional3D<T,3>, velocity.getBoundingBox(), velocity, S, envelopeWidth );
}

template<typename T>
std::auto_ptr<TensorField3D<T,3> > computeStrainRate(TensorField3D<T,3>& velocity)
{
    TensorField3D<T,6>* S = new TensorField3D<T,6>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeStrainRate(velocity, *S);
    return std::auto_ptr<TensorField3D<T,6> >(S);
}


/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkStrainRate(TensorField3D<T,3>& velocity, TensorField3D<T,6>& S)
{
    applyProcessingFunctional (
            new BoxBulkStrainRateFunctional3D<T,3>, velocity.getBoundingBox(), velocity, S );
}

template<typename T>
std::auto_ptr<TensorField3D<T,6> > computeBulkStrainRate(TensorField3D<T,3>& velocity)
{
    TensorField3D<T,6>* S = new TensorField3D<T,6>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeBulkStrainRate(velocity, *S);
    return std::auto_ptr<TensorField3D<T,6> >(S);
}


/* *************** TensorField - TensorField operations *************** */

template<typename T1, typename T2, int nDim>
void copy(TensorField3D<T1,nDim>& field, TensorField3D<T2,nDim>& convertedField)
{
    applyProcessingFunctional (
            new CopyConvertTensorFunctional3D<T1,T2,nDim>, field.getBoundingBox(), field, convertedField );
}

template<typename T1, typename T2, int nDim>
std::auto_ptr<TensorField3D<T2,nDim> > copyConvert(TensorField3D<T1,nDim>& field)
{
    TensorField3D<T2,nDim>* convertedField
        = new TensorField3D<T2,nDim>(field.getNx(), field.getNy(), field.getNz());
    plb::copy(field, *convertedField);
    return std::auto_ptr<TensorField3D<T2,nDim> >(convertedField);
}

template<typename T, int nDim>
void add(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B, TensorField3D<T,nDim>& result)
{
    std::vector<TensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_plus_B_functional3D<T,nDim>, A.getBoundingBox(), fields );
}

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > add(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    TensorField3D<T,nDim>* result = new TensorField3D<T,nDim>(A.getNx(), A.getNy(), A.getNz());
    add(A, B, *result);
    return std::auto_ptr<TensorField3D<T,nDim> >(result);
}


template<typename T, int nDim>
void subtract(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B, TensorField3D<T,nDim>& result)
{
    std::vector<TensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_minus_B_functional3D<T,nDim>, A.getBoundingBox(), fields );
}

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > subtract(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    TensorField3D<T,nDim>* result = new TensorField3D<T,nDim>(A.getNx(), A.getNy(), A.getNz());
    subtract(A, B, *result);
    return std::auto_ptr<TensorField3D<T,nDim> >(result);
}


template<typename T, int nDim>
void multiply(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B, TensorField3D<T,nDim>& result)
{
    std::vector<TensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_times_B_functional3D<T,nDim>, A.getBoundingBox(), fields );
}

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > multiply(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    TensorField3D<T,nDim>* result = new TensorField3D<T,nDim>(A.getNx(), A.getNy(), A.getNz());
    multiply(A, B, *result);
    return std::auto_ptr<TensorField3D<T,nDim> >(result);
}


template<typename T, int nDim>
void multiply(TensorField3D<T,nDim>& field, T scalar, TensorField3D<T,nDim>& result, Box3D domain)
{
    applyProcessingFunctional (
        new Tensor_A_times_alpha_functional3D<T,nDim>(scalar), domain, field, result );
}

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > multiply(TensorField3D<T,nDim>& field, T scalar)
{
    TensorField3D<T,nDim>* result = new TensorField3D<T,nDim>(field.getNx(), field.getNy(), field.getNz());
    multiply(field, scalar, *result);
    return std::auto_ptr<TensorField3D<T,nDim> >(result);
}


template<typename T, int nDim>
void multiply(T scalar, TensorField3D<T,nDim>& field, TensorField3D<T,nDim>& result, Box3D domain)
{
    applyProcessingFunctional (
        new Tensor_A_times_alpha_functional3D<T,nDim>(scalar), domain, field, result );
}

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > multiply(T scalar, TensorField3D<T,nDim>& field)
{
    TensorField3D<T,nDim>* result = new TensorField3D<T,nDim>(field.getNx(), field.getNy(), field.getNz());
    multiply(scalar, field, *result);
    return std::auto_ptr<TensorField3D<T,nDim> >(result);
}


template<typename T, int nDim>
void divide(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B, TensorField3D<T,nDim>& result)
{
    std::vector<TensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_dividedBy_B_functional3D<T,nDim>, A.getBoundingBox(), fields );
}

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > divide(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    TensorField3D<T,nDim>* result = new TensorField3D<T,nDim>(A.getNx(), A.getNy(), A.getNz());
    divide(A, B, *result);
    return std::auto_ptr<TensorField3D<T,nDim> >(result);
}


/* *************** TensorField - TensorField inplace operations *************** */

template<typename T, int nDim>
void addInPlace(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B) {
    applyProcessingFunctional (
            new Tensor_A_plus_B_inplace_functional3D<T,nDim>, A.getBoundingBox(), A, B );
}

template<typename T, int nDim>
void subtractInPlace(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B) {
    applyProcessingFunctional (
            new Tensor_A_minus_B_inplace_functional3D<T,nDim>, A.getBoundingBox(), A, B );
}

template<typename T, int nDim>
void multiplyInPlace(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B) {
    applyProcessingFunctional (
            new Tensor_A_times_B_inplace_functional3D<T,nDim>, A.getBoundingBox(), A, B );
}

template<typename T, int nDim>
void multiplyInPlace(TensorField3D<T,nDim>& A, T alpha)
{
    applyProcessingFunctional (
            new Tensor_A_times_alpha_inplace_functional3D<T,nDim>(alpha), A.getBoundingBox(), A);
}

template<typename T, int nDim>
void divideInPlace(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B) {
    applyProcessingFunctional (
            new Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>, A.getBoundingBox(), A, B );
}


/* ******************************************************************* */
/* *************** PART IV. Multi-block wrappers: Block-Lattice ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template<typename T, template<typename U> class Descriptor> 
T computeAverageDensity(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain) {
    BoxSumRhoBarFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return Descriptor<T>::fullRho( functional.getSumRhoBar() / (T) domain.nCells() );
}

template<typename T, template<typename U> class Descriptor> 
T computeAverageDensity(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computeAverageDensity(lattice, lattice.getBoundingBox());
}


template<typename T, template<typename U> class Descriptor> 
T computeAverageRhoBar(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain) {
    BoxSumRhoBarFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumRhoBar() / (T) domain.nCells();
}

template<typename T, template<typename U> class Descriptor> 
T computeAverageRhoBar(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computeAverageRhoBar(lattice, lattice.getBoundingBox());
}


template<typename T, template<typename U> class Descriptor> 
T computeAverageEnergy(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain) 
{
    BoxSumEnergyFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumEnergy() / (T) domain.nCells();;
}

template<typename T, template<typename U> class Descriptor> 
T computeAverageEnergy(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computeAverageEnergy(lattice, lattice.getBoundingBox());
}

template<typename T, template<typename U> class Descriptor, class BoolMask> 
plint count(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain, BoolMask boolMask)
{
    CountLatticeElementsFunctional3D<T,Descriptor,BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getCount();
}

template<typename T, template<typename U> class Descriptor, class BoolMask> 
plint count(MultiBlockLattice3D<T,Descriptor>& lattice, BoolMask boolMask)
{
    return count(lattice, lattice.getBoundingBox(), boolMask);
}

template<typename T, template<typename U> class Descriptor>
std::vector<T> densitySingleProbes (
        MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
        std::vector<Array<T,3> > const& positions )
{
    DensitySingleProbe3D<T,Descriptor> functional(positions);
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getDensities();
}

template<typename T, template<typename U> class Descriptor>
std::vector<T> densitySingleProbes (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        std::vector<Array<T,3> > const& positions )
{
    return densitySingleProbes(lattice, lattice.getBoundingBox(), positions);
}

template<typename T, template<typename U> class Descriptor>
std::vector<Array<T,3> > velocitySingleProbes (
        MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
        std::vector<Array<T,3> > const& positions )
{
    VelocitySingleProbe3D<T,Descriptor> functional(positions);;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getVelocities();
}

template<typename T, template<typename U> class Descriptor>
std::vector<Array<T,3> > velocitySingleProbes (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        std::vector<Array<T,3> > const& positions )
{
    return velocitySingleProbes(lattice, lattice.getBoundingBox(), positions);
}

template<typename T, template<typename U> class Descriptor>
std::vector<Array<T,3> > vorticitySingleProbes (
        MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
        std::vector<Array<T,3> > const& positions )
{
    VorticitySingleProbe3D<T,Descriptor> functional(positions);;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getVorticities();
}

template<typename T, template<typename U> class Descriptor>
std::vector<Array<T,3> > vorticitySingleProbes (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        std::vector<Array<T,3> > const& positions )
{
    return vorticitySingleProbes(lattice, lattice.getBoundingBox(), positions);
}

/* *************** Extract Sub-Lattice ******************************* */

template<typename T, template<typename U> class Descriptor>
void extractSubDomain( MultiBlockLattice3D<T,Descriptor>& lattice,
                       MultiBlockLattice3D<T,Descriptor>& extractedLattice,
                       Box3D domain)
{
    applyProcessingFunctional (
            new LatticeRegenerateFunctional3D<T,Descriptor>, domain, lattice, extractedLattice );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiBlockLattice3D<T,Descriptor> > extractSubDomain(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    std::auto_ptr<MultiBlockLattice3D<T,Descriptor> > extractedLattice =
        generateMultiBlockLattice<T,Descriptor>(lattice, domain);
    extractSubDomain(lattice, *extractedLattice, domain);
    return extractedLattice;
}


/* *************** Density ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeDensity(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& density, Box3D domain)
{
    applyProcessingFunctional (
            new BoxDensityFunctional3D<T,Descriptor>, domain, lattice, density );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeDensity(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > density =
        generateMultiScalarField<T>(lattice, domain);
    computeDensity(lattice, *density, domain);
    return density;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeDensity(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computeDensity(lattice, lattice.getBoundingBox());
}


/* *************** RhoBar ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeRhoBar(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& rhoBar, Box3D domain)
{
    applyProcessingFunctional (
            new BoxRhoBarFunctional3D<T,Descriptor>, domain, lattice, rhoBar );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeRhoBar(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > rhoBar =
        generateMultiScalarField<T>(lattice, domain);
    computeRhoBar(lattice, *rhoBar, domain);
    return rhoBar;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeRhoBar(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computeRhoBar(lattice, lattice.getBoundingBox());
}


template<typename T, template<typename U> class Descriptor>
void computeRhoBarJ( MultiBlockLattice3D<T,Descriptor>& lattice,
                     MultiScalarField3D<T>& rhoBar, MultiTensorField3D<T,3>& j, Box3D domain )
{
    std::vector<MultiBlock3D*> fields;
    fields.push_back(&lattice);
    fields.push_back(&rhoBar);
    fields.push_back(&j);
    applyProcessingFunctional (
            new BoxRhoBarJfunctional3D<T,Descriptor>, domain, fields );
}

/* *************** Kinetic Energy ************************************ */

template<typename T, template<typename U> class Descriptor>
void computeKineticEnergy(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& energy, Box3D domain)
{
    applyProcessingFunctional (
            new BoxKineticEnergyFunctional3D<T,Descriptor>, domain, lattice, energy );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeKineticEnergy(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > energy =
        generateMultiScalarField<T>(lattice, domain);
    computeKineticEnergy(lattice, *energy, domain);
    return energy;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeKineticEnergy(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computeKineticEnergy(lattice, lattice.getBoundingBox());
}


/* *************** Packed RhoBar J *********************************** */

template<typename T, template<typename U> class Descriptor>
void computePackedRhoBarJ(MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<T>& rhoBarJ, Box3D domain)
{
    applyProcessingFunctional (
            new PackedRhoBarJfunctional3D<T,Descriptor>, domain, lattice, rhoBarJ );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiNTensorField3D<T> > computePackedRhoBarJ(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    std::auto_ptr<MultiNTensorField3D<T> > rhoBarJ =
        generateMultiNTensorField<T>(lattice, domain);
    computePackedRhoBarJ(lattice, *rhoBarJ, domain);
    return rhoBarJ;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computePackedRhoBarJ(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computePackedRhoBarJ(lattice, lattice.getBoundingBox());
}


template<typename T>
void computeDensityFromRhoBarJ (
        MultiNTensorField3D<T>& rhoBarJ,
        MultiScalarField3D<T>& density, Box3D domain)
{
    applyProcessingFunctional (
            new DensityFromRhoBarJfunctional3D<T>, domain, density, rhoBarJ );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeDensityFromRhoBarJ (
        MultiNTensorField3D<T>& rhoBarJ, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > density =
        generateMultiScalarField<T>(rhoBarJ, domain);
    computeDensityFromRhoBarJ(rhoBarJ, *density, domain);
    return density;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeDensityFromRhoBarJ(MultiNTensorField3D<T>& rhoBarJ) {
    return computeDensityFromRhoBarJ(rhoBarJ, rhoBarJ.getBoundingBox());
}


template<typename T>
void computeVelocityFromRhoBarJ (
        MultiNTensorField3D<T>& rhoBarJ,
        MultiTensorField3D<T,3>& velocity, Box3D domain, bool velIsJ )
{
    std::vector<MultiBlock3D*> args;
    args.push_back(&velocity);
    args.push_back(&rhoBarJ);
    applyProcessingFunctional (
            new VelocityFromRhoBarJfunctional3D<T>(velIsJ), domain, args );
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeVelocityFromRhoBarJ (
        MultiNTensorField3D<T>& rhoBarJ, Box3D domain, bool velIsJ )
{
    std::auto_ptr<MultiTensorField3D<T,3> > velocity =
        generateMultiTensorField<T,3>(rhoBarJ, domain);
    computeVelocityFromRhoBarJ(rhoBarJ, *velocity, domain, velIsJ);
    return velocity;
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeVelocityFromRhoBarJ(MultiNTensorField3D<T>& rhoBarJ, bool velIsJ)
{
    return computeVelocityFromRhoBarJ(rhoBarJ, rhoBarJ.getBoundingBox(), velIsJ);
}

/* *************** Velocity Norm ************************************* */

template<typename T, template<typename U> class Descriptor>
void computeVelocityNorm(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& velocityNorm, Box3D domain)
{
    applyProcessingFunctional (
            new BoxVelocityNormFunctional3D<T,Descriptor>, domain, lattice, velocityNorm );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeVelocityNorm(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > velocityNorm =
        generateMultiScalarField<T>(lattice, domain);
    computeVelocityNorm(lattice, *velocityNorm, domain);
    return velocityNorm;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeVelocityNorm(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computeVelocityNorm(lattice, lattice.getBoundingBox());
}


/* *************** Velocity Component ******************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocityComponent(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& velocityComponent,
                              Box3D domain, plint iComponent)
{
    applyProcessingFunctional (
            new BoxVelocityComponentFunctional3D<T,Descriptor>(iComponent), domain, lattice, velocityComponent );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeVelocityComponent(MultiBlockLattice3D<T,Descriptor>& lattice,
                                                               Box3D domain, plint iComponent)
{
    std::auto_ptr<MultiScalarField3D<T> > velocityComponent =
        generateMultiScalarField<T>(lattice, domain);
    computeVelocityComponent(lattice, *velocityComponent, domain, iComponent);
    return velocityComponent;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeVelocityComponent(MultiBlockLattice3D<T,Descriptor>& lattice, plint iComponent)
{
    return computeVelocityComponent(lattice, lattice.getBoundingBox(), iComponent);
}


/* *************** Velocity ****************************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocity(MultiBlockLattice3D<T,Descriptor>& lattice,
                     MultiTensorField3D<T,Descriptor<T>::d>& velocity, Box3D domain)
{
    applyProcessingFunctional (
            new BoxVelocityFunctional3D<T,Descriptor>, domain, lattice, velocity );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::d> > computeVelocity(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::d> > velocity
        = generateMultiTensorField<T,Descriptor<T>::d>(lattice, domain);
    computeVelocity(lattice, *velocity, domain);
    return velocity;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::d> >
    computeVelocity(MultiBlockLattice3D<T,Descriptor>& lattice)
{
    return computeVelocity(lattice, lattice.getBoundingBox());
}

/* *************** Temperature ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeTemperature(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& temperature, Box3D domain)
{
    applyProcessingFunctional (
        new BoxTemperatureFunctional3D<T,Descriptor>, domain, lattice, temperature );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeTemperature(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > temperature =
        generateMultiScalarField<T>(lattice, domain);
    computeTemperature(lattice, *temperature, domain);
    return temperature;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeTemperature(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computeTemperature(lattice, lattice.getBoundingBox());
}


/* *************** Pi Neq ********************************* */

template<typename T, template<typename U> class Descriptor>
void computePiNeq( MultiBlockLattice3D<T,Descriptor>& lattice,
                              MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq,
                              Box3D domain )
{
    applyProcessingFunctional (
            new BoxPiNeqFunctional3D<T,Descriptor>, domain, lattice, PiNeq );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
    computePiNeq(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> > PiNeq
        = generateMultiTensorField<T,SymmetricTensor<T,Descriptor>::n>(lattice, domain);
    computePiNeq(lattice, *PiNeq, domain);
    return PiNeq;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
    computePiNeq(MultiBlockLattice3D<T,Descriptor>& lattice)
{
    return computePiNeq(lattice, lattice.getBoundingBox());
}


/* *************** Shear Stress ********************************* */

template<typename T, template<typename U> class Descriptor>
void computeShearStress( MultiBlockLattice3D<T,Descriptor>& lattice,
                              MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq,
                              Box3D domain )
{
    applyProcessingFunctional (
            new BoxShearStressFunctional3D<T,Descriptor>, domain, lattice, PiNeq );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
    computeShearStress(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> > PiNeq
        = generateMultiTensorField<T,SymmetricTensor<T,Descriptor>::n>(lattice, domain);
    computeShearStress(lattice, *PiNeq, domain);
    return PiNeq;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
    computeShearStress(MultiBlockLattice3D<T,Descriptor>& lattice)
{
    return computeShearStress(lattice, lattice.getBoundingBox());
}


/* *************** Strain Rate from Stress *************************** */

template<typename T, template<typename U> class Descriptor>
void computeStrainRateFromStress( MultiBlockLattice3D<T,Descriptor>& lattice,
                                  MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n>& S,
                                  Box3D domain )
{
    applyProcessingFunctional (
            new BoxStrainRateFromStressFunctional3D<T,Descriptor>, domain, lattice, S );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> > S
        = generateMultiTensorField<T,SymmetricTensor<T,Descriptor>::n>(lattice, domain);
    computeStrainRateFromStress(lattice, *S, domain);
    return S;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice3D<T,Descriptor>& lattice)
{
    return computeStrainRateFromStress(lattice, lattice.getBoundingBox());
}


/* *************** Population **************************************** */

template<typename T, template<typename U> class Descriptor>
void computePopulation(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& population,
                       Box3D domain, plint iPop)
{
    applyProcessingFunctional (
            new BoxPopulationFunctional3D<T,Descriptor>(iPop), domain, lattice, population );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computePopulation(MultiBlockLattice3D<T,Descriptor>& lattice,
                                                        Box3D domain, plint iPop)
{
    std::auto_ptr<MultiScalarField3D<T> > population = generateMultiScalarField<T>(lattice, domain);
    computePopulation(lattice, *population, domain, iPop);
    return population;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computePopulation(MultiBlockLattice3D<T,Descriptor>& lattice, plint iPop)
{
    return computePopulation(lattice, lattice.getBoundingBox(), iPop);
}


template<typename T, template<typename U> class Descriptor>
void computeEquilibrium(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& equilibrium,
                       Box3D domain, plint iPop)
{
    applyProcessingFunctional (
            new BoxEquilibriumFunctional3D<T,Descriptor>(iPop), domain, lattice, equilibrium );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeEquilibrium(MultiBlockLattice3D<T,Descriptor>& lattice,
                                                        Box3D domain, plint iPop)
{
    std::auto_ptr<MultiScalarField3D<T> > equilibrium = generateMultiScalarField<T>(lattice, domain);
    computeEquilibrium(lattice, *equilibrium, domain, iPop);
    return equilibrium;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeEquilibrium(MultiBlockLattice3D<T,Descriptor>& lattice, plint iPop)
{
    return computeEquilibrium(lattice, lattice.getBoundingBox(), iPop);
}

template<typename T, template<typename U> class Descriptor>
void computeEquilibrium(MultiBlockLattice3D<T,Descriptor>& lattice, MultiTensorField3D<T,Descriptor<T>::q>& equilibrium,
                        Box3D domain)
{
    applyProcessingFunctional (
        new BoxAllEquilibriumFunctional3D<T,Descriptor>(), domain, lattice, equilibrium );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::q> > computeEquilibrium(MultiBlockLattice3D<T,Descriptor>& lattice,
                                                                          Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::q> > equilibrium = generateMultiTensorField<T,Descriptor<T>::q>(lattice, domain);
    computeEquilibrium(lattice, *equilibrium, domain);
    return equilibrium;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::q> > computeEquilibrium(MultiBlockLattice3D<T,Descriptor>& lattice)
{
    return computeEquilibrium(lattice, lattice.getBoundingBox());
}

template<typename T, template<typename U> class Descriptor>
void computeAllPopulations(MultiBlockLattice3D<T,Descriptor>& lattice, 
                           MultiTensorField3D<T,Descriptor<T>::q>& populations,
                           Box3D domain)
{
    applyProcessingFunctional (
            new BoxAllPopulationsFunctional3D<T,Descriptor>(), domain, lattice, populations );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::q> > computeAllPopulations(MultiBlockLattice3D<T,Descriptor>& lattice,
                                                            Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::q> > populations = 
            generateMultiTensorField<T,Descriptor<T>::q>(lattice, domain);
    computeAllPopulations(lattice, *populations, domain);
    return populations;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::q> > computeAllPopulations(MultiBlockLattice3D<T,Descriptor>& lattice)
{
    return computeAllPopulations(lattice, lattice.getBoundingBox());
}

template<typename T, template<typename U> class Descriptor>
void computeAllPopulationsFromTensorField(MultiBlockLattice3D<T,Descriptor>& lattice, 
                           MultiTensorField3D<T,Descriptor<T>::q>& populations,
                           Box3D domain)
{
    applyProcessingFunctional (
            new BoxAllPopulationsToLatticeFunctional3D<T,Descriptor>(), domain, lattice, populations );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::q> > computeAllPopulationsFromTensorField(MultiBlockLattice3D<T,Descriptor>& lattice,
                                                            Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::q> > populations = 
            generateMultiTensorField<T,Descriptor<T>::q>(lattice, domain);
    computeAllPopulationsFromTensorField(lattice, *populations, domain);
    return populations;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::q> > computeAllPopulationsFromTensorField(MultiBlockLattice3D<T,Descriptor>& lattice)
{
    return computeAllPopulationsFromTensorField(lattice, lattice.getBoundingBox());
}

template<typename T, template<typename U> class Descriptor>
void copyPopulations(MultiBlockLattice3D<T,Descriptor>& latticeFrom, MultiBlockLattice3D<T,Descriptor>& latticeTo,
                     Box3D domain)
{
    applyProcessingFunctional (
            new CopyPopulationsFunctional3D<T,Descriptor>(), domain, latticeFrom, latticeTo );
}

template<typename T, template<typename U> class Descriptor>
void copyAll(MultiBlockLattice3D<T,Descriptor>& latticeFrom,
             MultiBlockLattice3D<T,Descriptor>& latticeTo, Box3D domain)
{
    applyProcessingFunctional (
            new LatticeCopyAllFunctional3D<T,Descriptor>(), domain, latticeFrom, latticeTo );
}

template<typename T, template<typename U> class Descriptor>
void copyRegenerate(MultiBlockLattice3D<T,Descriptor>& latticeFrom,
                    MultiBlockLattice3D<T,Descriptor>& latticeTo, Box3D domain)
{
    applyProcessingFunctional (
            new LatticeRegenerateFunctional3D<T,Descriptor>(), domain, latticeFrom, latticeTo );
}


/* ******************************************************************* */
/* *************** PART V. Multi-block wrappers: Scalar-Field ******** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */
template<typename T>
T computeSum(MultiScalarField3D<T>& scalarField, Box3D domain) {
    BoxScalarSumFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar();
}

template<typename T>
T computeSum(MultiScalarField3D<T>& scalarField) {
    return computeSum(scalarField, scalarField.getBoundingBox());
}

template<typename T>
T computeBoundedSum(MultiScalarField3D<T>& scalarField, Box3D domain) {
    BoundedBoxScalarSumFunctional3D<T> functional;
    plint envelopeWidth=1;
    applyProcessingFunctional(functional, domain, scalarField, envelopeWidth);
    return functional.getSumScalar();
}

template<typename T>
T computeBoundedSum(MultiScalarField3D<T>& scalarField) {
    return computeBoundedSum(scalarField, scalarField.getBoundingBox());
}

template<typename T>
T computeAverage(MultiScalarField3D<T>& scalarField, Box3D domain) {
    BoxScalarSumFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar() / (T) domain.nCells();
}

template<typename T>
T computeAverage(MultiScalarField3D<T>& scalarField) {
    return computeAverage(scalarField, scalarField.getBoundingBox());
}


template<typename T>
T computeAverage(MultiScalarField3D<T>& scalarField, MultiScalarField3D<int>& mask, int flag, Box3D domain)
{
    MaskedBoxScalarAverageFunctional3D<T> functional(flag);
    applyProcessingFunctional(functional, domain, scalarField, mask);
    return functional.getAverageScalar();
}

template<typename T>
T computeAverage(MultiScalarField3D<T>& scalarField, MultiScalarField3D<int>& mask, int flag) {
    return computeAverage(scalarField, mask, flag, scalarField.getBoundingBox());
}


template<typename T>
T computeMin(MultiScalarField3D<T>& scalarField, Box3D domain) {
    BoxScalarMinFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getMinScalar();
}

template<typename T>
T computeMin(MultiScalarField3D<T>& scalarField) {
    return computeMin(scalarField, scalarField.getBoundingBox());
}


template<typename T>
T computeMax(MultiScalarField3D<T>& scalarField, Box3D domain) {
    BoxScalarMaxFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getMaxScalar();
}

template<typename T>
T computeMax(MultiScalarField3D<T>& scalarField) {
    return computeMax(scalarField, scalarField.getBoundingBox());
}


template<typename T>
T computeBoundedAverage(MultiScalarField3D<T>& scalarField, Box3D domain) {
    BoundedBoxScalarSumFunctional3D<T> functional;
    plint envelopeWidth=1;
    applyProcessingFunctional(functional, domain, scalarField, envelopeWidth);
    return functional.getSumScalar() /
             (T) ( (domain.getNx()-1)*(domain.getNy()-1)*(domain.getNz()-1) );
}

template<typename T>
T computeBoundedAverage(MultiScalarField3D<T>& scalarField) {
    return computeBoundedAverage(scalarField, scalarField.getBoundingBox());
}


template<typename T, class BoolMask> 
plint count(MultiScalarField3D<T>& field, Box3D domain, BoolMask boolMask)
{
    CountScalarElementsFunctional3D<T,BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

template<typename T, class BoolMask> 
plint count(MultiScalarField3D<T>& field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

/* *************** Extract Sub-ScalarField *************************** */

template<typename T>
void extractSubDomain( MultiScalarField3D<T>& field,
                       MultiScalarField3D<T>& extractedField,
                       Box3D domain)
{
    applyProcessingFunctional (
            new ExtractScalarSubDomainFunctional3D<T>, domain, field, extractedField );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > extractSubDomain(MultiScalarField3D<T>& field, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > extractedField =
        generateMultiScalarField<T>(field, domain);
    extractSubDomain(field, *extractedField, domain);
    return extractedField;
}

/* *************** MultiScalarField - Scalar operations *************** */


template<typename T>
void lessThan(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<int>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_lt_alpha_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<int> > lessThan(MultiScalarField3D<T>& field, T scalar, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<int> > result =
        generateMultiScalarField<int>(field, domain);
    lessThan(field, scalar, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<int> > lessThan(MultiScalarField3D<T>& field, T scalar)
{
    return lessThan(field, scalar, field.getBoundingBox());
}


template<typename T>
void greaterThan(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<int>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_gt_alpha_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<int> > greaterThan(MultiScalarField3D<T>& field, T scalar, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<int> > result =
        generateMultiScalarField<int>(field, domain);
    greaterThan(field, scalar, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<int> > greaterThan(MultiScalarField3D<T>& field, T scalar)
{
    return greaterThan(field, scalar, field.getBoundingBox());
}

template<typename T>
void add(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_plus_alpha_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T>& field, T scalar, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(field, domain);
    add(field, scalar, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T>& field, T scalar)
{
    return add(field, scalar, field.getBoundingBox());
}


template<typename T>
void add(T scalar, MultiScalarField3D<T>& field, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_plus_alpha_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(T scalar, MultiScalarField3D<T>& field, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(field, domain);
    add(scalar, field, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(T scalar, MultiScalarField3D<T>& field)
{
    return add(scalar, field, field.getBoundingBox());
}


template<typename T>
void subtract(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_minus_alpha_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T>& field, T scalar, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(field, domain);
    subtract(field, scalar, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T>& field, T scalar)
{
    return subtract(field, scalar, field.getBoundingBox());
}


template<typename T>
void subtract(T scalar, MultiScalarField3D<T>& field, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new Alpha_minus_A_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(T scalar, MultiScalarField3D<T>& field, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(field, domain);
    subtract(scalar, field, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(T scalar, MultiScalarField3D<T>& field)
{
    return subtract(scalar, field, field.getBoundingBox());
}


template<typename T>
void multiply(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_times_alpha_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T>& field, T scalar, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(field, domain);
    multiply(field, scalar, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T>& field, T scalar)
{
    return multiply(field, scalar, field.getBoundingBox());
}


template<typename T>
void multiply(T scalar, MultiScalarField3D<T>& field, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_times_alpha_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(T scalar, MultiScalarField3D<T>& field, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(field, domain);
    multiply(scalar, field, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(T scalar, MultiScalarField3D<T>& field)
{
    return multiply(scalar, field, field.getBoundingBox());
}


template<typename T>
void divide(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_dividedBy_alpha_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T>& field, T scalar, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(field, domain);
    divide(field, scalar, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T>& field, T scalar)
{
    return divide(field, scalar, field.getBoundingBox());
}


template<typename T>
void divide(T scalar, MultiScalarField3D<T>& field, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new Alpha_dividedBy_A_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(T scalar, MultiScalarField3D<T>& field, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(field, domain);
    divide(scalar, field, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(T scalar, MultiScalarField3D<T>& field)
{
    return divide(scalar, field, field.getBoundingBox());
}


/* *************** MultiScalarField - Scalar inplace operations *************** */

template<typename T>
void addInPlace(MultiScalarField3D<T>& field, T scalar, Box3D domain) {
    applyProcessingFunctional (
            new A_plus_alpha_inplace_functional3D<T>(scalar), domain, field);
}

template<typename T>
void addInPlace(MultiScalarField3D<T>& field, T scalar) {
    addInPlace(field, scalar, field.getBoundingBox());
}


template<typename T>
void subtractInPlace(MultiScalarField3D<T>& field, T scalar, Box3D domain) {
    applyProcessingFunctional (
            new A_minus_alpha_inplace_functional3D<T>(scalar), domain, field);
}

template<typename T>
void subtractInPlace(MultiScalarField3D<T>& field, T scalar) {
    subtractInPlace(field, scalar, field.getBoundingBox());
}


template<typename T>
void multiplyInPlace(MultiScalarField3D<T>& field, T scalar, Box3D domain) {
    applyProcessingFunctional (
            new A_times_alpha_inplace_functional3D<T>(scalar), domain, field);
}

template<typename T>
void multiplyInPlace(MultiScalarField3D<T>& field, T scalar) {
    multiplyInPlace(field, scalar, field.getBoundingBox());
}


template<typename T>
void divideInPlace(MultiScalarField3D<T>& field, T scalar, Box3D domain) {
    applyProcessingFunctional (
            new A_dividedBy_alpha_inplace_functional3D<T>(scalar), domain, field);
}

template<typename T>
void divideInPlace(MultiScalarField3D<T>& field, T scalar) {
    divideInPlace(field, scalar, field.getBoundingBox());
}


/* *************** MultiScalarField - MultiScalarField operations *************** */

template<typename T1, typename T2>
void copy(MultiScalarField3D<T1>& field, MultiScalarField3D<T2>& convertedField, Box3D domain)
{
    applyProcessingFunctional (
            new CopyConvertScalarFunctional3D<T1,T2>, domain, field, convertedField );
}

template<typename T1, typename T2>
std::auto_ptr<MultiScalarField3D<T2> > copyConvert(MultiScalarField3D<T1>& field, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T2> > convertedField =
        generateMultiScalarField<T2>(field, domain);
    plb::copy(field, *convertedField, domain);
    return convertedField;
}

template<typename T1, typename T2>
std::auto_ptr<MultiScalarField3D<T2> > copyConvert(MultiScalarField3D<T1>& field)
{
    return copyConvert<T1,T2>(field, field.getBoundingBox());
}


template<typename T>
void lessThan(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<int>& result, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_lt_B_functional3D<T>, domain, fields );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<int> > lessThan(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<int> > result =
        generateIntersectMultiScalarField<int>(A,B, domain);
    lessThan(A, B, *result, domain);
    return result;
}


template<typename T>
std::auto_ptr<MultiScalarField3D<int> > lessThan(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B)
{
    return lessThan(A, B, A.getBoundingBox());
}


template<typename T>
void greaterThan(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<int>& result, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_gt_B_functional3D<T>, domain, fields );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<int> > greaterThan(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<int> > result =
        generateIntersectMultiScalarField<int>(A,B, domain);
    greaterThan(A, B, *result, domain);
    return result;
}


template<typename T>
std::auto_ptr<MultiScalarField3D<int> > greaterThan(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B)
{
    return greaterThan(A, B, A.getBoundingBox());
}



template<typename T>
void add(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<T>& result, Box3D domain)
{
    std::vector<MultiScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_plus_B_functional3D<T>, domain, fields );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateIntersectMultiScalarField<T>(A,B, domain);
    add(A, B, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B)
{
    return add(A, B, A.getBoundingBox());
}


template<typename T>
void subtract(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<T>& result, Box3D domain)
{
    std::vector<MultiScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_minus_B_functional3D<T>, domain, fields );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateIntersectMultiScalarField<T>(A,B, domain);
    subtract(A, B, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B)
{
    return subtract(A, B, A.getBoundingBox());
}


template<typename T>
void multiply(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<T>& result, Box3D domain)
{
    std::vector<MultiScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_times_B_functional3D<T>, domain, fields );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateIntersectMultiScalarField<T>(A,B, domain);
    multiply(A, B, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B)
{
    return multiply(A, B, A.getBoundingBox());
}

template<typename T>
void divide(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<T>& result, Box3D domain)
{
    std::vector<MultiScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_dividedBy_B_functional3D<T>, domain, fields );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateIntersectMultiScalarField<T>(A,B, domain);
    divide(A, B, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B)
{
    return divide(A, B, A.getBoundingBox());
}


/* *************** ScalarField - ScalarField inplace operations *************** */

template<typename T>
void addInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain) {
    applyProcessingFunctional (
            new A_plus_B_inplace_functional3D<T>, domain, A, B );
}

template<typename T>
void addInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B) {
    addInPlace(A, B, A.getBoundingBox());
}


template<typename T>
void subtractInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain) {
    applyProcessingFunctional (
            new A_minus_B_inplace_functional3D<T>, domain, A, B );
}

template<typename T>
void subtractInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B) {
    subtractInPlace(A, B, A.getBoundingBox());
}


template<typename T>
void multiplyInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain) {
    applyProcessingFunctional (
            new A_times_B_inplace_functional3D<T>, domain, A, B );
}

template<typename T>
void multiplyInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B) {
    multiplyInPlace(A,B, A.getBoundingBox());
}


template<typename T>
void divideInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain) {
    applyProcessingFunctional (
            new A_dividedBy_B_inplace_functional3D<T>, domain, A, B );
}

template<typename T>
void divideInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B) {
    divideInPlace(A, B, A.getBoundingBox());
}

/* *************** MultiScalarField operations *************** */

template<typename T>
void computeSqrt(MultiScalarField3D<T>& A, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeScalarSqrtFunctional3D<T>, domain, A, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSqrt(MultiScalarField3D<T>& A, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(A, domain);
    computeSqrt(A, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSqrt(MultiScalarField3D<T>& A)
{
    return computeSqrt(A, A.getBoundingBox());
}


template<typename T>
void computeAbsoluteValue(MultiScalarField3D<T>& A, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
        new ComputeAbsoluteValueFunctional3D<T>, domain, A, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeAbsoluteValue(MultiScalarField3D<T>& A, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
    generateMultiScalarField<T>(A, domain);
    computeAbsoluteValue(A, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeAbsoluteValue(MultiScalarField3D<T>& A)
{
    return computeAbsoluteValue(A, A.getBoundingBox());
}


template<typename T>
void uniformlyBoundScalarField(MultiScalarField3D<T>& data, Box3D domain, T bound) {
    applyProcessingFunctional (
            new UniformlyBoundScalarField3D<T>(bound), domain, data );
}

template<typename T>
void uniformlyBoundScalarField(MultiScalarField3D<T>& data, T bound) {
    uniformlyBoundScalarField(data, data.getBoundingBox(), bound);
}


/* ******************************************************************* */
/* *************** PART VI. Multi-block wrappers: Tensor-field ******* */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template<typename T, plint nDim, class BoolMask> 
plint count(MultiTensorField3D<T,nDim>& field, Box3D domain, BoolMask boolMask)
{
    CountTensorElementsFunctional3D<T,nDim,BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

template<typename T, plint nDim, class BoolMask> 
plint count(MultiTensorField3D<T,nDim>& field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

/* *************** Extract Sub-TensorField *************************** */

template<typename T, int nDim>
void extractSubDomain( MultiTensorField3D<T,nDim>& field,
                       MultiTensorField3D<T,nDim>& extractedField,
                       Box3D domain)
{
    applyProcessingFunctional (
            new ExtractTensorSubDomainFunctional3D<T,nDim>, domain, field, extractedField );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > extractSubDomain(MultiTensorField3D<T,nDim>& field, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,nDim> > extractedField =
        generateMultiTensorField<T,nDim>(field, domain);
    extractSubDomain(field, *extractedField, domain);
    return extractedField;
}


/* *************** Component (scalar-field) out of a tensor-field ****** */

template<typename T, int nDim>
void extractComponent(MultiTensorField3D<T,nDim>& tensorField, MultiScalarField3D<T>& component, Box3D domain, int iComponent)
{
    applyProcessingFunctional (
            new ExtractTensorComponentFunctional3D<T,nDim>(iComponent), domain, component, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > extractComponent(MultiTensorField3D<T,nDim>& tensorField, Box3D domain, int iComponent)
{
    std::auto_ptr<MultiScalarField3D<T> > component =
        generateMultiScalarField<T>(tensorField, domain);
    extractComponent(tensorField, *component, domain, iComponent);
    return component;
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > extractComponent(MultiTensorField3D<T,nDim>& tensorField, int iComponent)
{
    return extractComponent(tensorField, tensorField.getBoundingBox(), iComponent);
}


/* *************** Vector-norm of each cell in the field *************** */

template<typename T, int nDim>
void computeNorm(MultiTensorField3D<T,nDim>& tensorField, MultiScalarField3D<T>& norm, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeNormFunctional3D<T,nDim>, domain, norm, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > computeNorm(MultiTensorField3D<T,nDim>& tensorField, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > norm = generateMultiScalarField<T>(tensorField, domain);
    computeNorm(tensorField, *norm, domain);
    return norm;
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > computeNorm(MultiTensorField3D<T,nDim>& tensorField)
{
    return computeNorm(tensorField, tensorField.getBoundingBox());
}

/* *************** Sqrt operation on each component of each cell *************** */

template<typename T, int nDim>
void computeSqrt(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeTensorSqrtFunctional3D<T, nDim>, domain, A, result);
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > computeSqrt(MultiTensorField3D<T,nDim>& A, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,nDim> > result =
        generateMultiTensorField<T,nDim>(A, domain);
    computeSqrt(A, *result, domain);
    return result;
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > computeSqrt(MultiTensorField3D<T,nDim>& A)
{
    return computeSqrt(A, A.getBoundingBox());
}


/* *************** Squared vector-norm of each cell in the field ******** */

template<typename T, int nDim>
void computeNormSqr(MultiTensorField3D<T,nDim>& tensorField, MultiScalarField3D<T>& normSqr, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeNormSqrFunctional3D<T,nDim>, domain, normSqr, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > computeNormSqr(MultiTensorField3D<T,nDim>& tensorField, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > normSqr =
        generateMultiScalarField<T>(tensorField, domain);
    computeNormSqr(tensorField, *normSqr, domain);
    return normSqr;
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > computeNormSqr(MultiTensorField3D<T,nDim>& tensorField)
{
    return computeNormSqr(tensorField, tensorField.getBoundingBox());
}


/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template<typename T>
void computeSymmetricTensorNorm(MultiTensorField3D<T,6>& tensorField, MultiScalarField3D<T>& norm, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorNormFunctional3D<T>, domain, norm, tensorField );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorNorm(MultiTensorField3D<T,6>& tensorField, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > norm =
        generateMultiScalarField<T>(tensorField, domain);
    computeSymmetricTensorNorm(tensorField, *norm, domain);
    return norm;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorNorm(MultiTensorField3D<T,6>& tensorField)
{
    return computeSymmetricTensorNorm(tensorField, tensorField.getBoundingBox());
}


/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template<typename T>
void computeSymmetricTensorNormSqr(MultiTensorField3D<T,6>& tensorField, MultiScalarField3D<T>& normSqr, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorNormSqrFunctional3D<T>, domain, normSqr, tensorField );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorNormSqr(MultiTensorField3D<T,6>& tensorField, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > normSqr =
        generateMultiScalarField<T>(tensorField, domain);
    computeSymmetricTensorNormSqr(tensorField, *normSqr, domain);
    return normSqr;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorNormSqr(MultiTensorField3D<T,6>& tensorField)
{
    return computeSymmetricTensorNormSqr(tensorField, tensorField.getBoundingBox());
}


/* *************** Trace of each symmetric tensor of a field ************ */

template<typename T>
void computeSymmetricTensorTrace(MultiTensorField3D<T,6>& tensorField, MultiScalarField3D<T>& trace, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorTraceFunctional3D<T>, domain, trace, tensorField );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorTrace(MultiTensorField3D<T,6>& tensorField, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > trace =
        generateMultiScalarField<T>(tensorField, domain);
    computeSymmetricTensorTrace(tensorField, *trace, domain);
    return trace;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorTrace(MultiTensorField3D<T,6>& tensorField)
{
    return computeSymmetricTensorTrace(tensorField, tensorField.getBoundingBox());
}


/* *************** Gradient from scalar field *********************** */

template<typename T>
void computeGradient(MultiScalarField3D<T>& phi, MultiTensorField3D<T,3>& gradient, Box3D domain)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxGradientFunctional3D<T>, domain, phi, gradient, envelopeWidth );
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeGradient(MultiScalarField3D<T>& phi, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,3> > gradient =
        generateMultiTensorField<T,3>(phi, domain);
    computeGradient(phi, *gradient, domain);
    return gradient;
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeGradient(MultiScalarField3D<T>& phi)
{
    return computeGradient(phi, phi.getBoundingBox());
}


/* *************** Gradient, without boundary treatment, from scalar field  */

template<typename T>
void computeBulkGradient(MultiScalarField3D<T>& phi, MultiTensorField3D<T,3>& gradient, Box3D domain)
{
    applyProcessingFunctional (
            new BoxBulkGradientFunctional3D<T>, domain, phi, gradient );
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeBulkGradient(MultiScalarField3D<T>& phi, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,3> > gradient =
        generateMultiTensorField<T,3>(phi, domain);
    computeBulkGradient(phi, *gradient, domain);
    return gradient;
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeBulkGradient(MultiScalarField3D<T>& phi)
{
    return computeBulkGradient(phi, phi.getBoundingBox());
}



/* *************** Vorticity from Velocity field *********************** */

template<typename T>
void computeVorticity(MultiTensorField3D<T,3>& velocity, MultiTensorField3D<T,3>& vorticity, Box3D domain)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxVorticityFunctional3D<T,3>, domain, velocity, vorticity, envelopeWidth );
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeVorticity(MultiTensorField3D<T,3>& velocity, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,3> > vorticity =
        generateMultiTensorField<T,3>(velocity, domain);
    computeVorticity(velocity, *vorticity, domain);
    return vorticity;
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeVorticity(MultiTensorField3D<T,3>& velocity)
{
    return computeVorticity(velocity, velocity.getBoundingBox());
}


/* *************** Vorticity, without boundary treatment, from Velocity field  */

template<typename T>
void computeBulkVorticity(MultiTensorField3D<T,3>& velocity, MultiTensorField3D<T,3>& vorticity, Box3D domain)
{
    applyProcessingFunctional (
            new BoxBulkVorticityFunctional3D<T,3>, domain, velocity, vorticity );
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeBulkVorticity(MultiTensorField3D<T,3>& velocity, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,3> > vorticity =
        generateMultiTensorField<T,3>(velocity, domain);
    computeBulkVorticity(velocity, *vorticity, domain);
    return vorticity;
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeBulkVorticity(MultiTensorField3D<T,3>& velocity)
{
    return computeBulkVorticity(velocity, velocity.getBoundingBox());
}



/* *************** Divergence, without boundary treatment, from Velocity field  */

template<typename T>
void computeBulkDivergence(MultiTensorField3D<T,3>& velocity, MultiScalarField3D<T>& divergence, Box3D domain)
{
    applyProcessingFunctional (
            new BoxBulkDivergenceFunctional3D<T,3>, domain, divergence, velocity );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeBulkDivergence(MultiTensorField3D<T,3>& velocity, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > divergence =
        generateMultiScalarField<T>(velocity, domain);
    computeBulkDivergence(velocity, *divergence, domain);
    return divergence;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeBulkDivergence(MultiTensorField3D<T,3>& velocity)
{
    return computeBulkDivergence(velocity, velocity.getBoundingBox());
}



/* *************** Strain rate from Velocity field ********************* */

template<typename T>
void computeStrainRate(MultiTensorField3D<T,3>& velocity, MultiTensorField3D<T,6>& S, Box3D domain)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxStrainRateFunctional3D<T,3>, domain, velocity, S, envelopeWidth );
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,6> > computeStrainRate(MultiTensorField3D<T,3>& velocity, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,6> > S =
        generateMultiTensorField<T,6>(velocity, domain);
    computeStrainRate(velocity, *S, domain);
    return S;
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,6> > computeStrainRate(MultiTensorField3D<T,3>& velocity)
{
    return computeStrainRate(velocity, velocity.getBoundingBox());
}


/* *************** Str. rate, without boundary treatment, from Velocity field  */

template<typename T>
void computeBulkStrainRate(MultiTensorField3D<T,3>& velocity, MultiTensorField3D<T,6>& S, Box3D domain)
{
    applyProcessingFunctional (
            new BoxBulkStrainRateFunctional3D<T,3>, domain, velocity, S );
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,6> > computeBulkStrainRate(MultiTensorField3D<T,3>& velocity, Box3D domain)
{
     std::auto_ptr<MultiTensorField3D<T,6> > S =
         generateMultiTensorField<T,6>(velocity, domain);
    computeBulkStrainRate(velocity, *S, domain);
    return S;
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,6> > computeBulkStrainRate(MultiTensorField3D<T,3>& velocity)
{
    return computeBulkStrainRate(velocity, velocity.getBoundingBox());
}

/* *************** Q-criterion from vorticity and strain rate fields  */

template<typename T>
void computeQcriterion(MultiTensorField3D<T,3>& vorticity, MultiTensorField3D<T,6>& S, MultiScalarField3D<T>& qCriterion, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&vorticity);
    fields.push_back(&S);
    fields.push_back(&qCriterion);
    
    applyProcessingFunctional (
        new BoxQcriterionFunctional3D<T>, domain, fields  );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeQcriterion(MultiTensorField3D<T,3>& vorticity, MultiTensorField3D<T,6>& S, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > qCriterion = generateMultiScalarField<T>(vorticity, domain);
    computeQcriterion(vorticity, S, *qCriterion, domain);
    return qCriterion;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeQcriterion(MultiTensorField3D<T,3>& vorticity, MultiTensorField3D<T,6>& S)
{
    return computeQcriterion(vorticity, S, vorticity.getBoundingBox());
}


/* *************** MultiTensorField - MultiTensorField operations *************** */

template<typename T1, typename T2, int nDim>
void copy(MultiTensorField3D<T1,nDim>& field, MultiTensorField3D<T2,nDim>& convertedField, Box3D domain)
{
    applyProcessingFunctional (
            new CopyConvertTensorFunctional3D<T1,T2,nDim>, domain, field, convertedField );
}

template<typename T1, typename T2, int nDim>
std::auto_ptr<MultiTensorField3D<T2,nDim> > copyConvert(MultiTensorField3D<T1,nDim>& field, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T2,nDim> > convertedField =
        generateMultiTensorField<T2,nDim>(field, domain);
    plb::copy<T1,T2,nDim>(field, *convertedField, domain);
    return convertedField;
}

template<typename T1, typename T2, int nDim>
std::auto_ptr<MultiTensorField3D<T2,nDim> > copyConvert(MultiTensorField3D<T1,nDim>& field)
{
    return copyConvert<T1,T2,nDim>(field, field.getBoundingBox());
}

template<typename T, int nDim>
void add(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    std::vector<MultiTensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_plus_B_functional3D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > add(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,nDim> > result =
        generateIntersectMultiTensorField<T,nDim>(A,B, domain);
    add(A, B, *result, domain);
    return result;
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > add(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B)
{
    return add(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void subtract(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    std::vector<MultiTensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_minus_B_functional3D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > subtract(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,nDim> > result =
        generateIntersectMultiTensorField<T,nDim>(A,B, domain);
    subtract(A, B, *result, domain);
    return result;
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > subtract(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B)
{
    return subtract(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void multiply(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    std::vector<MultiTensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_times_B_functional3D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > multiply(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,nDim> > result =
        generateIntersectMultiTensorField<T,nDim>(A,B, domain);
    multiply(A, B, *result, domain);
    return result;
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > multiply(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B)
{
    return multiply(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void multiply(MultiTensorField3D<T,nDim>& field, T scalar, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    applyProcessingFunctional (
            new Tensor_A_times_alpha_functional3D<T,nDim>(scalar), domain, field, result );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > multiply(MultiTensorField3D<T,nDim>& field, T scalar, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,nDim> > result =
        generateMultiTensorField<T,nDim>(field, domain);
    multiply(field, scalar, *result, domain);
    return result;
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > multiply(MultiTensorField3D<T,nDim>& field, T scalar)
{
    return multiply(field, scalar, field.getBoundingBox());
}


template<typename T, int nDim>
void multiply(T scalar, MultiTensorField3D<T,nDim>& field, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    applyProcessingFunctional (
            new Tensor_A_times_alpha_functional3D<T,nDim>(scalar), domain, field, result );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > multiply(T scalar, MultiTensorField3D<T,nDim>& field, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,nDim> > result =
        generateMultiTensorField<T,nDim>(field, domain);
    multiply(scalar, field, *result, domain);
    return result;
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > multiply(T scalar, MultiTensorField3D<T,nDim>& field)
{
    return multiply(scalar, field, field.getBoundingBox());
}


template<typename T, int nDim>
void divide(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    std::vector<MultiTensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
    new Tensor_A_dividedBy_B_functional3D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > divide(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,nDim> > result =
        generateIntersectMultiTensorField<T,nDim>(A,B, domain);
    divide(A, B, *result, domain);
    return result;
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > divide(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B)
{
    return divide(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void normalize(MultiTensorField3D<T,nDim>& data, MultiTensorField3D<T,nDim>& result, Box3D domain, Precision precision)
{
    applyProcessingFunctional (
            new Normalize_Tensor_functional3D<T,nDim>(precision), domain, data, result );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > normalize(MultiTensorField3D<T,nDim>& data, Box3D domain, Precision precision)
{
    std::auto_ptr<MultiTensorField3D<T,nDim> > result =
        generateMultiTensorField<T,nDim>(data, domain);
    normalize(data, *result, domain, precision);
    return result;
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > normalize(MultiTensorField3D<T,nDim>& data, Precision precision)
{
    return normalize(data, data.getBoundingBox(), precision);
}


// tensor product of a vector with itself
template<typename T, int nDim>
void symmetricTensorProduct(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& result, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&result);
    applyProcessingFunctional (
            new TensorProduct_A_A_functional3D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > symmetricTensorProduct(MultiTensorField3D<T,nDim>& A, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,SymmetricTensorImpl<T,nDim>::n> > result =
        generateMultiTensorField<T,nDim>(A,domain);
    symmetricTensorProduct(A, *result, domain);
    return result;
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > symmetricTensorProduct(MultiTensorField3D<T,nDim>& A)
{
    return symmetricTensorProduct(A, A.getBoundingBox());
}


/* *************** MultiScalarField - MultiTensorField operations *************** */

template<typename T, int nDim>
void multiply(MultiScalarField3D<T>& A, MultiTensorField3D<T,nDim>& B, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Scalar_A_times_Tensor_B_functional3D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > multiply(MultiScalarField3D<T>& A, MultiTensorField3D<T,nDim>& B, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,nDim> > result =
        generateMultiTensorField<T,nDim>(B, domain);
    multiply(A, B, *result, domain);
    return result;
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > multiply(MultiScalarField3D<T>& A, MultiTensorField3D<T,nDim>& B)
{
    return multiply(A, B, A.getBoundingBox());
}



template<typename T, int nDim>
void fullIndexContractionOfSymmetricTensors(MultiTensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& A,
                                            MultiTensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& B,
                                            MultiScalarField3D<T>& result, Box3D domain)
{
    std::vector<MultiBlock3D* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > fullIndexContractionOfSymmetricTensors(
    MultiTensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& A,
    MultiTensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& B, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(B, domain);
    fullIndexContractionOfSymmetricTensors(A, B, *result, domain);
    return result;
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > fullIndexContractionOfSymmetricTensors(
    MultiTensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& A, MultiTensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& B)
{
    return fullIndexContractionOfSymmetricTensors(A, B, A.getBoundingBox());
}



/* *************** TensorField - TensorField inplace operations *************** */

template<typename T, int nDim>
void addInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain) {
    applyProcessingFunctional (
            new Tensor_A_plus_B_inplace_functional3D<T,nDim>, domain, A, B );
}

template<typename T, int nDim>
void addInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B) {
    addInPlace(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void subtractInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain) {
    applyProcessingFunctional (
            new Tensor_A_minus_B_inplace_functional3D<T,nDim>, domain, A, B );
}

template<typename T, int nDim>
void subtractInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B) {
    subtractInPlace(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain) {
    applyProcessingFunctional (
            new Tensor_A_times_B_inplace_functional3D<T,nDim>, domain, A, B );
}

template<typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B) {
    multiplyInPlace(A,B, A.getBoundingBox());
}


template<typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T,nDim>& A, T alpha, Box3D domain) {
    applyProcessingFunctional (
            new Tensor_A_times_alpha_inplace_functional3D<T,nDim>(alpha), domain, A);
}

template<typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T,nDim>& A, T alpha)
{
    multiplyInPlace(A,alpha, A.getBoundingBox());
}


template<typename T, int nDim>
void divideInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain) {
    applyProcessingFunctional (
            new Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>, domain, A, B );
}

template<typename T, int nDim>
void divideInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B) {
    divideInPlace(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void normalizeInPlace(MultiTensorField3D<T,nDim>& data, Box3D domain, Precision precision) {
    applyProcessingFunctional (
            new Normalize_Tensor_inplace_functional3D<T,nDim>(precision), domain, data );
}

template<typename T, int nDim>
void normalizeInPlace(MultiTensorField3D<T,nDim>& data, Precision precision) {
    normalizeInPlace(data, data.getBoundingBox(), precision);
}


/* *************** LBMsmoothen3D ******************************************* */

template<typename T, template<typename U> class Descriptor>
void lbmSmoothen(MultiScalarField3D<T>& data, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional(new LBMsmoothen3D<T,Descriptor>, domain, data, result);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > lbmSmoothen(MultiScalarField3D<T>& data, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(data, domain);
    lbmSmoothen<T,Descriptor>(data, *result, domain);
    return result;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > lbmSmoothen(MultiScalarField3D<T>& data)
{
    return lbmSmoothen<T,Descriptor>(data, data.getBoundingBox());
}


/* *************** Smoothen3D ******************************************* */

template<typename T>
void smoothen(MultiScalarField3D<T>& data, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional(new Smoothen3D<T>, domain, data, result);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > smoothen(MultiScalarField3D<T>& data, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(data, domain);
    smoothen<T>(data, *result, domain);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > smoothen(MultiScalarField3D<T>& data)
{
    return smoothen<T>(data, data.getBoundingBox());
}

/* *************** LBMsmoothenTensor3D ******************************************* */

template<typename T, int nDim, template<typename U> class Descriptor>
void lbmSmoothenTensor(MultiTensorField3D<T,nDim>& data, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    applyProcessingFunctional(new LBMsmoothenTensor3D<T,nDim,Descriptor>, domain, data, result);
}

template<typename T, int nDim, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,nDim> > lbmSmoothenTensor(MultiTensorField3D<T,nDim>& data, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,nDim> > result =
        generateMultiTensorField<T,nDim>(data, domain);
    lbmSmoothenTensor<T,nDim,Descriptor>(data, *result, domain);
    return result;
}

template<typename T, int nDim, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,nDim> > lbmSmoothenTensor(MultiTensorField3D<T,nDim>& data)
{
    return lbmSmoothenTensor<T,nDim,Descriptor>(data, data.getBoundingBox());
}

/* *************** SmoothenTensor3D ******************************************* */

template<typename T, int nDim>
void smoothenTensor(MultiTensorField3D<T,nDim>& data, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    applyProcessingFunctional(new SmoothenTensor3D<T,nDim>, domain, data, result);
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > smoothenTensor(MultiTensorField3D<T,nDim>& data, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,nDim> > result =
        generateMultiTensorField<T,nDim>(data, domain);
    smoothenTensor<T,nDim>(data, *result, domain);
    return result;
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > smoothenTensor(MultiTensorField3D<T,nDim>& data)
{
    return smoothenTensor<T,nDim>(data, data.getBoundingBox());
}


/* *************** MollifyScalar3D ******************************************* */

template<typename T>
void mollifyScalar(T l, plint d, Box3D globalDomain, int exclusionFlag,
        MultiScalarField3D<T>& data, MultiScalarField3D<int>& flag, MultiScalarField3D<T>& result, Box3D domain)
{
    std::vector<MultiBlock3D*> args;
    args.push_back(&data);
    args.push_back(&flag);
    args.push_back(&result);
    applyProcessingFunctional(new MollifyScalar3D<T>(l, d, globalDomain, exclusionFlag), domain, args);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > mollifyScalar(T l, plint d, Box3D globalDomain, int exclusionFlag,
        MultiScalarField3D<T>& data, MultiScalarField3D<int>& flag, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(data, domain);

    std::vector<MultiBlock3D*> args;
    args.push_back(&data);
    args.push_back(&flag);
    args.push_back(result.get());
    applyProcessingFunctional(new MollifyScalar3D<T>(l, d, globalDomain, exclusionFlag), domain, args);

    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > mollifyScalar(T l, plint d, Box3D globalDomain, int exclusionFlag,
        MultiScalarField3D<T>& data, MultiScalarField3D<int>& flag)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(data, data.getBoundingBox());

    std::vector<MultiBlock3D*> args;
    args.push_back(&data);
    args.push_back(&flag);
    args.push_back(result.get());
    applyProcessingFunctional(new MollifyScalar3D<T>(l, d, globalDomain, exclusionFlag), data.getBoundingBox(), args);

    return result;
}


/* *************** MollifyTensor3D ******************************************* */

template<typename T, int nDim>
void mollifyTensor(T l, plint d, Box3D globalDomain, int exclusionFlag,
        MultiTensorField3D<T,nDim>& data, MultiScalarField3D<int>& flag, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    std::vector<MultiBlock3D*> args;
    args.push_back(&data);
    args.push_back(&flag);
    args.push_back(&result);
    applyProcessingFunctional(new MollifyTensor3D<T,nDim>(l, d, globalDomain, exclusionFlag), domain, args);
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > mollifyTensor(T l, plint d, Box3D globalDomain, int exclusionFlag,
        MultiTensorField3D<T,nDim>& data, MultiScalarField3D<int>& flag, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,nDim> > result =
        generateMultiTensorField<T,nDim>(data, domain);

    std::vector<MultiBlock3D*> args;
    args.push_back(&data);
    args.push_back(&flag);
    args.push_back(result.get());
    applyProcessingFunctional(new MollifyTensor3D<T,nDim>(l, d, globalDomain, exclusionFlag), domain, args);

    return result;
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > mollifyTensor(T l, plint d, Box3D globalDomain, int exclusionFlag,
        MultiTensorField3D<T,nDim>& data, MultiScalarField3D<int>& flag)
{
    std::auto_ptr<MultiTensorField3D<T,nDim> > result =
        generateMultiTensorField<T,nDim>(data, data.getBoundingBox());

    std::vector<MultiBlock3D*> args;
    args.push_back(&data);
    args.push_back(&flag);
    args.push_back(result.get());
    applyProcessingFunctional(new MollifyTensor3D<T,nDim>(l, d, globalDomain, exclusionFlag), data.getBoundingBox(), args);

    return result;
}


/* *************** LBMcomputeGradient3D ******************************************* */

template<typename T, template<typename U> class Descriptor>
void lbmComputeGradient(MultiScalarField3D<T>& scalarField, MultiTensorField3D<T,3>& gradient, Box3D domain)
{
    applyProcessingFunctional(new LBMcomputeGradient3D<T,Descriptor>, domain, scalarField, gradient);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,3> > lbmComputeGradient(MultiScalarField3D<T>& scalarField, Box3D domain)
{
    std::auto_ptr<MultiTensorField3D<T,3> > gradient =
        generateMultiTensorField<T,3>(scalarField, domain);
    lbmComputeGradient<T,Descriptor>(scalarField, *gradient, domain);
    return gradient;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,3> > lbmComputeGradient(MultiScalarField3D<T>& scalarField) {
    return lbmComputeGradient<T,Descriptor>(scalarField, scalarField.getBoundingBox());
}

/* *************** LBMcomputeDivergence3D ******************************************* */

template<typename T, template<typename U> class Descriptor>
void lbmComputeDivergence(MultiScalarField3D<T>& divergence, MultiTensorField3D<T,3>& vectorField, Box3D domain)
{
    applyProcessingFunctional(new LBMcomputeDivergence3D<T,Descriptor>, domain, divergence, vectorField);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > lbmComputeDivergence(MultiTensorField3D<T,3>& vectorField, Box3D domain)
{
    std::auto_ptr<MultiScalarField3D<T> > divergence =
        generateMultiScalarField<T>(vectorField, domain);
    lbmComputeDivergence<T,Descriptor>(*divergence, vectorField, domain);
    return divergence;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > lbmComputeDivergence(MultiTensorField3D<T,3>& vectorField) {
    return lbmComputeDivergence<T,Descriptor>(vectorField, vectorField.getBoundingBox());
}

}  // namespace plb

#endif  // DATA_ANALYSIS_WRAPPER_3D_HH

