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

#ifdef COMPILE_3D

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- template instantiation.
 */
#include "plbWrapper/lattice/latticeAnalysisWrapper3D.h"
#include "plbWrapper/lattice/latticeAnalysisWrapper3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

template void computeDensity<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& density, Box3D domain );

template void maskedComputeDensity<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& density, MultiNTensorField3D<int>& mask, Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* pypal_computeDensity (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* maskedPypal_computeDensity (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain );

template void computeKineticEnergy<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& energy, Box3D domain );

template void maskedComputeKineticEnergy<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& energy, MultiNTensorField3D<int>& mask, Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* pypal_computeKineticEnergy (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* maskedPypal_computeKineticEnergy (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain );


template void computeVelocityNorm<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& velocityNorm, Box3D domain );

template void maskedComputeVelocityNorm<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& velocityNorm, MultiNTensorField3D<int>& mask, Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* pypal_computeVelocityNorm (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* maskedPypal_computeVelocityNorm (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain );


template void computeVelocityComponent<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& velocityComponent, Box3D domain, plint iComponent );

template void maskedComputeVelocityComponent<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& velocityComponent, MultiNTensorField3D<int>& mask, Box3D domain, plint iComponent );

template
MultiNTensorField3D<FLOAT_T>* pypal_computeVelocityComponent (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        Box3D domain, plint iComponent );

template
MultiNTensorField3D<FLOAT_T>* maskedPypal_computeVelocityComponent (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain, plint iComponent );

template void computeVelocity<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& velocity, Box3D domain );

template void maskedComputeVelocity<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& velocity, MultiNTensorField3D<int>& mask, Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* pypal_computeVelocity (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* maskedPypal_computeVelocity (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
void computePiNeq (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& PiNeq, Box3D domain );

template
void maskedComputePiNeq (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& PiNeq, MultiNTensorField3D<int>& mask, Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* pypal_computePiNeq (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* maskedPypal_computePiNeq (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain );

template
void computeShearStress (
    MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
    MultiNTensorField3D<FLOAT_T>& ShearStress, Box3D domain );

template
void maskedComputeShearStress (
    MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
    MultiNTensorField3D<FLOAT_T>& ShearStress, MultiNTensorField3D<int>& mask, Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* pypal_computeShearStress (
    MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
    Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* maskedPypal_computeShearStress (
    MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice, MultiNTensorField3D<int>& mask,
    Box3D domain );

template
void computeStrainRateFromStress (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& S, Box3D domain );

template
void maskedComputeStrainRateFromStress (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& S, MultiNTensorField3D<int>& mask, Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* pypal_computeStrainRateFromStress (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* maskedPypal_computeStrainRateFromStress (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain );

template void computePopulation<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& population, Box3D domain, plint iPop );

template void maskedComputePopulation<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& population, MultiNTensorField3D<int>& mask, Box3D domain, plint iPop );

template
MultiNTensorField3D<FLOAT_T>* pypal_computePopulation (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        Box3D domain, plint iPop );

template
MultiNTensorField3D<FLOAT_T>* maskedPypal_computePopulation (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain, plint iPop );


template void computePopulations<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& populations, Box3D domain );

template void maskedComputePopulations<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& populations, MultiNTensorField3D<int>& mask, Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* pypal_computePopulations (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* maskedPypal_computePopulations (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain );


template
void compute_UPO_Rhs<FLOAT_T,descriptors::DESCRIPTOR_3D>(
                     MultiNTensorField3D<FLOAT_T>& lattice,
                     MultiNTensorField3D<FLOAT_T>& result, Box3D domain, FLOAT_T omega);

template
void masked_compute_UPO_Rhs<FLOAT_T,descriptors::DESCRIPTOR_3D>(
                            MultiNTensorField3D<FLOAT_T>& lattice,
                            MultiNTensorField3D<FLOAT_T>& result,
                            MultiNTensorField3D<int>& mask, Box3D domain, FLOAT_T omega);

template
MultiNTensorField3D<FLOAT_T>* compute_UPO_Rhs<FLOAT_T,descriptors::DESCRIPTOR_3D> (
                           MultiNTensorField3D<FLOAT_T>& lattice, Box3D domain, FLOAT_T omega );

template
MultiNTensorField3D<FLOAT_T>* masked_compute_UPO_Rhs<FLOAT_T,descriptors::DESCRIPTOR_3D> (
                                 MultiNTensorField3D<FLOAT_T>& lattice,
                                 MultiNTensorField3D<int>& mask, Box3D domain, FLOAT_T omega );



template
void UPO_ApplyJ<FLOAT_T,descriptors::DESCRIPTOR_3D>(
                 MultiNTensorField3D<FLOAT_T>& f,
                 MultiNTensorField3D<FLOAT_T>& g,
                 MultiNTensorField3D<FLOAT_T>& h, Box3D domain, FLOAT_T omega);

template
void masked_UPO_ApplyJ<FLOAT_T,descriptors::DESCRIPTOR_3D>(
                        MultiNTensorField3D<FLOAT_T>& f,
                        MultiNTensorField3D<FLOAT_T>& g,
                        MultiNTensorField3D<FLOAT_T>& h,
                        MultiNTensorField3D<int>& mask, Box3D domain, FLOAT_T omega);

template
MultiNTensorField3D<FLOAT_T>* UPO_ApplyJ<FLOAT_T,descriptors::DESCRIPTOR_3D> (
                           MultiNTensorField3D<FLOAT_T>& f,
                           MultiNTensorField3D<FLOAT_T>& g,
                           Box3D domain, FLOAT_T omega );

template
MultiNTensorField3D<FLOAT_T>* masked_UPO_ApplyJ<FLOAT_T,descriptors::DESCRIPTOR_3D> (
                                 MultiNTensorField3D<FLOAT_T>& f,
                                 MultiNTensorField3D<FLOAT_T>& g,
                                 MultiNTensorField3D<int>& mask, Box3D domain, FLOAT_T omega );



template
void compute_UPO_EnergyDerivative<FLOAT_T,descriptors::DESCRIPTOR_3D> (
                 MultiNTensorField3D<FLOAT_T>& lattice,
                 MultiNTensorField3D<FLOAT_T>& result, Box3D domain );

template
void masked_compute_UPO_EnergyDerivative<FLOAT_T,descriptors::DESCRIPTOR_3D> (
                        MultiNTensorField3D<FLOAT_T>& lattice,
                        MultiNTensorField3D<FLOAT_T>& result,
                        MultiNTensorField3D<int>& mask, Box3D domain);

template
MultiNTensorField3D<FLOAT_T>* compute_UPO_EnergyDerivative<FLOAT_T,descriptors::DESCRIPTOR_3D> (
                           MultiNTensorField3D<FLOAT_T>& lattice, Box3D domain );

template
MultiNTensorField3D<FLOAT_T>* masked_compute_UPO_EnergyDerivative<FLOAT_T,descriptors::DESCRIPTOR_3D> (
                                 MultiNTensorField3D<FLOAT_T>& lattice,
                                 MultiNTensorField3D<int>& mask, Box3D domain );

}  // namespace plb

#endif  // COMPILE_3D
