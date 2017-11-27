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

#ifdef COMPILE_2D

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- template instantiation.
 */
#include "plbWrapper/lattice/latticeAnalysisWrapper2D.h"
#include "plbWrapper/lattice/latticeAnalysisWrapper2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {

template void computeDensity<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& density, Box2D domain );

template void maskedComputeDensity<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& density, MultiNTensorField2D<int>& mask, Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* pypal_computeDensity (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* maskedPypal_computeDensity (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain );

template void computeKineticEnergy<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& energy, Box2D domain );

template void maskedComputeKineticEnergy<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& energy, MultiNTensorField2D<int>& mask, Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* pypal_computeKineticEnergy (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* maskedPypal_computeKineticEnergy (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain );


template void computeVelocityNorm<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& velocityNorm, Box2D domain );

template void maskedComputeVelocityNorm<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& velocityNorm, MultiNTensorField2D<int>& mask, Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* pypal_computeVelocityNorm (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* maskedPypal_computeVelocityNorm (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain );


template void computeVelocityComponent<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& velocityComponent, Box2D domain, plint iComponent );

template void maskedComputeVelocityComponent<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& velocityComponent, MultiNTensorField2D<int>& mask, Box2D domain, plint iComponent );

template
MultiNTensorField2D<FLOAT_T>* pypal_computeVelocityComponent (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        Box2D domain, plint iComponent );

template
MultiNTensorField2D<FLOAT_T>* maskedPypal_computeVelocityComponent (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain, plint iComponent );

template void computeVelocity<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& velocity, Box2D domain );

template void maskedComputeVelocity<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& velocity, MultiNTensorField2D<int>& mask, Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* pypal_computeVelocity (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* maskedPypal_computeVelocity (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain );


template
void computePiNeq (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& PiNeq, Box2D domain );

template
void maskedComputePiNeq (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& PiNeq, MultiNTensorField2D<int>& mask, Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* pypal_computePiNeq (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* maskedPypal_computePiNeq (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain );

template
void computeShearStress (
    MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
    MultiNTensorField2D<FLOAT_T>& ShearStress, Box2D domain );

template
void maskedComputeShearStress (
    MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
    MultiNTensorField2D<FLOAT_T>& ShearStress, MultiNTensorField2D<int>& mask, Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* pypal_computeShearStress (
    MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
    Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* maskedPypal_computeShearStress (
    MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice, MultiNTensorField2D<int>& mask,
    Box2D domain );

template
void computeStrainRateFromStress (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& S, Box2D domain );

template
void maskedComputeStrainRateFromStress (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& S, MultiNTensorField2D<int>& mask, Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* pypal_computeStrainRateFromStress (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* maskedPypal_computeStrainRateFromStress (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain );

template void computePopulation<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& population, Box2D domain, plint iPop );

template void maskedComputePopulation<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& population, MultiNTensorField2D<int>& mask, Box2D domain, plint iPop );

template
MultiNTensorField2D<FLOAT_T>* pypal_computePopulation (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        Box2D domain, plint iPop );

template
MultiNTensorField2D<FLOAT_T>* maskedPypal_computePopulation (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain, plint iPop );


template void computePopulations<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& populations, Box2D domain );

template void maskedComputePopulations<FLOAT_T, descriptors::DESCRIPTOR_2D> (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        MultiNTensorField2D<FLOAT_T>& populations, MultiNTensorField2D<int>& mask, Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* pypal_computePopulations (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice,
        Box2D domain );

template
MultiNTensorField2D<FLOAT_T>* maskedPypal_computePopulations (
        MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain );

}  // namespace plb

#endif  // COMPILE_2D
