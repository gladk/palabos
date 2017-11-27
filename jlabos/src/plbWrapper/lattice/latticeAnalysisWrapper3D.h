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
 * Helper functions for domain initialization -- header file.
 */
#ifndef SWIG_LATTICE_ANALYSIS_WRAPPER_3D_H
#define SWIG_LATTICE_ANALYSIS_WRAPPER_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include <memory>


namespace plb {

template<typename T, template<typename U> class Descriptor>
void computeDensity(MultiBlockLattice3D<T,Descriptor>& lattice,
                    MultiNTensorField3D<T>& density, Box3D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputeDensity(MultiBlockLattice3D<T,Descriptor>& lattice,
                          MultiNTensorField3D<T>& density, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computeDensity (
        MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computeDensity (
              MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask, Box3D domain );

template<typename T, template<typename U> class Descriptor>
void computeKineticEnergy(MultiBlockLattice3D<T,Descriptor>& lattice,
                          MultiNTensorField3D<T>& energy, Box3D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputeKineticEnergy(MultiBlockLattice3D<T,Descriptor>& lattice,
                                MultiNTensorField3D<T>& energy, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computeKineticEnergy (
                            MultiBlockLattice3D<T,Descriptor>& lattice,
                            Box3D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computeKineticEnergy (
                                  MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
                                  Box3D domain );

template<typename T, template<typename U> class Descriptor>
void computeVelocityNorm(MultiBlockLattice3D<T,Descriptor>& lattice,
                         MultiNTensorField3D<T>& velocityNorm, Box3D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputeVelocityNorm(MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiNTensorField3D<T>& velocityNorm, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computeVelocityNorm (
                            MultiBlockLattice3D<T,Descriptor>& lattice,
                            Box3D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computeVelocityNorm (
                                  MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
                                  Box3D domain );

template<typename T, template<typename U> class Descriptor>
void computeVelocityComponent(MultiBlockLattice3D<T,Descriptor>& lattice,
                              MultiNTensorField3D<T>& velocityComponent,
                              Box3D domain, plint iComponent);

template<typename T, template<typename U> class Descriptor>
void maskedComputeVelocityComponent(MultiBlockLattice3D<T,Descriptor>& lattice,
                                    MultiNTensorField3D<T>& velocityComponent, MultiNTensorField3D<int>& mask,
                                    Box3D domain, plint iComponent);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computeVelocityComponent (
                            MultiBlockLattice3D<T,Descriptor>& lattice,
                            Box3D domain, plint iComponent );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computeVelocityComponent (
                                  MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
                                  Box3D domain, plint iComponent );

template<typename T, template<typename U> class Descriptor>
void computeVelocity(MultiBlockLattice3D<T,Descriptor>& lattice,
                     MultiNTensorField3D<T>& velocity, Box3D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputeVelocity(MultiBlockLattice3D<T,Descriptor>& lattice,
                           MultiNTensorField3D<T>& velocity, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computeVelocity (
                            MultiBlockLattice3D<T,Descriptor>& lattice,
                            Box3D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computeVelocity (
                                  MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
                                  Box3D domain );

template<typename T, template<typename U> class Descriptor>
void computePiNeq(MultiBlockLattice3D<T,Descriptor>& lattice,
                             MultiNTensorField3D<T>& PiNeq, Box3D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputePiNeq(MultiBlockLattice3D<T,Descriptor>& lattice,
                                   MultiNTensorField3D<T>& PiNeq, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computePiNeq (
                            MultiBlockLattice3D<T,Descriptor>& lattice,
                            Box3D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computePiNeq (
                                  MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
                                  Box3D domain );


template<typename T, template<typename U> class Descriptor>
void computeShearStress(MultiBlockLattice3D<T,Descriptor>& lattice,
                  MultiNTensorField3D<T>& ShearStress, Box3D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputeShearStress(MultiBlockLattice3D<T,Descriptor>& lattice,
                        MultiNTensorField3D<T>& ShearStress, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computeShearStress (
    MultiBlockLattice3D<T,Descriptor>& lattice,
    Box3D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computeShearStress (
    MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
    Box3D domain );

template<typename T, template<typename U> class Descriptor>
void computeStrainRateFromStress(MultiBlockLattice3D<T,Descriptor>& lattice,
                                 MultiNTensorField3D<T>& S, Box3D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputeStrainRateFromStress(MultiBlockLattice3D<T,Descriptor>& lattice,
                                       MultiNTensorField3D<T>& S, MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computeStrainRateFromStress (
                           MultiBlockLattice3D<T,Descriptor>& lattice,
                           Box3D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computeStrainRateFromStress (
                                 MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
                                 Box3D domain );

template<typename T, template<typename U> class Descriptor>
void computePopulation(MultiBlockLattice3D<T,Descriptor>& lattice,
                       MultiNTensorField3D<T>& population,
                       Box3D domain, plint iPop);

template<typename T, template<typename U> class Descriptor>
void maskedComputePopulation(MultiBlockLattice3D<T,Descriptor>& lattice,
                             MultiNTensorField3D<T>& population, MultiNTensorField3D<int>& mask,
                             Box3D domain, plint iPop);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computePopulation (
                            MultiBlockLattice3D<T,Descriptor>& lattice,
                            Box3D domain, plint iPop );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computePopulation (
                                  MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
                                  Box3D domain, plint iPop );


template<typename T, template<typename U> class Descriptor>
void computePopulations(MultiBlockLattice3D<T,Descriptor>& lattice,
                        MultiNTensorField3D<T>& populations,
                        Box3D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputePopulations(MultiBlockLattice3D<T,Descriptor>& lattice,
                              MultiNTensorField3D<T>& populations, MultiNTensorField3D<int>& mask,
                              Box3D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computePopulations (
                            MultiBlockLattice3D<T,Descriptor>& lattice,
                            Box3D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computePopulations (
                                  MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
                                  Box3D domain );


template<typename T, template<typename U> class Descriptor>
void compute_UPO_Rhs(MultiNTensorField3D<T>& lattice,
                     MultiNTensorField3D<T>& result, Box3D domain, T omega);

template<typename T, template<typename U> class Descriptor>
void masked_compute_UPO_Rhs(MultiNTensorField3D<T>& lattice,
                            MultiNTensorField3D<T>& result,
                            MultiNTensorField3D<int>& mask, Box3D domain, T omega);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* compute_UPO_Rhs (
                           MultiNTensorField3D<T>& lattice, Box3D domain, T omega );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* masked_compute_UPO_Rhs (
                                 MultiNTensorField3D<T>& lattice,
                                 MultiNTensorField3D<int>& mask, Box3D domain, T omega );



template<typename T, template<typename U> class Descriptor>
void UPO_ApplyJ( MultiNTensorField3D<T>& f,
                 MultiNTensorField3D<T>& g,
                 MultiNTensorField3D<T>& h, Box3D domain, T omega);

template<typename T, template<typename U> class Descriptor>
void masked_UPO_ApplyJ( MultiNTensorField3D<T>& f,
                        MultiNTensorField3D<T>& g,
                        MultiNTensorField3D<T>& h,
                        MultiNTensorField3D<int>& mask, Box3D domain, T omega);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* UPO_ApplyJ (
                           MultiNTensorField3D<T>& f,
                           MultiNTensorField3D<T>& g,
                           Box3D domain, T omega );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* masked_UPO_ApplyJ (
                                 MultiNTensorField3D<T>& f,
                                 MultiNTensorField3D<T>& g,
                                 MultiNTensorField3D<int>& mask, Box3D domain, T omega );



template<typename T, template<typename U> class Descriptor>
void compute_UPO_EnergyDerivative (
                 MultiNTensorField3D<T>& lattice,
                 MultiNTensorField3D<T>& result, Box3D domain );

template<typename T, template<typename U> class Descriptor>
void masked_compute_UPO_EnergyDerivative (
                        MultiNTensorField3D<T>& lattice,
                        MultiNTensorField3D<T>& result,
                        MultiNTensorField3D<int>& mask, Box3D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* compute_UPO_EnergyDerivative (
                           MultiNTensorField3D<T>& lattice, Box3D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* masked_compute_UPO_EnergyDerivative (
                                 MultiNTensorField3D<T>& lattice,
                                 MultiNTensorField3D<int>& mask, Box3D domain );


}  // namespace plb

#endif  // SWIG_LATTICE_ANALYSIS_WRAPPER_3D_H
