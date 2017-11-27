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
#ifndef SWIG_LATTICE_ANALYSIS_WRAPPER_2D_H
#define SWIG_LATTICE_ANALYSIS_WRAPPER_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
#include <memory>


namespace plb {

template<typename T, template<typename U> class Descriptor>
void computeDensity(MultiBlockLattice2D<T,Descriptor>& lattice,
                    MultiNTensorField2D<T>& density, Box2D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputeDensity(MultiBlockLattice2D<T,Descriptor>& lattice,
                          MultiNTensorField2D<T>& density, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computeDensity (
        MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computeDensity (
              MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask, Box2D domain );

template<typename T, template<typename U> class Descriptor>
void computeKineticEnergy(MultiBlockLattice2D<T,Descriptor>& lattice,
                          MultiNTensorField2D<T>& energy, Box2D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputeKineticEnergy(MultiBlockLattice2D<T,Descriptor>& lattice,
                                MultiNTensorField2D<T>& energy, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computeKineticEnergy (
                            MultiBlockLattice2D<T,Descriptor>& lattice,
                            Box2D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computeKineticEnergy (
                                  MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
                                  Box2D domain );

template<typename T, template<typename U> class Descriptor>
void computeVelocityNorm(MultiBlockLattice2D<T,Descriptor>& lattice,
                         MultiNTensorField2D<T>& velocityNorm, Box2D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputeVelocityNorm(MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiNTensorField2D<T>& velocityNorm, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computeVelocityNorm (
                            MultiBlockLattice2D<T,Descriptor>& lattice,
                            Box2D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computeVelocityNorm (
                                  MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
                                  Box2D domain );

template<typename T, template<typename U> class Descriptor>
void computeVelocityComponent(MultiBlockLattice2D<T,Descriptor>& lattice,
                              MultiNTensorField2D<T>& velocityComponent,
                              Box2D domain, plint iComponent);

template<typename T, template<typename U> class Descriptor>
void maskedComputeVelocityComponent(MultiBlockLattice2D<T,Descriptor>& lattice,
                                    MultiNTensorField2D<T>& velocityComponent, MultiNTensorField2D<int>& mask,
                                    Box2D domain, plint iComponent);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computeVelocityComponent (
                            MultiBlockLattice2D<T,Descriptor>& lattice,
                            Box2D domain, plint iComponent );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computeVelocityComponent (
                                  MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
                                  Box2D domain, plint iComponent );

template<typename T, template<typename U> class Descriptor>
void computeVelocity(MultiBlockLattice2D<T,Descriptor>& lattice,
                     MultiNTensorField2D<T>& velocity, Box2D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputeVelocity(MultiBlockLattice2D<T,Descriptor>& lattice,
                           MultiNTensorField2D<T>& velocity, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computeVelocity (
                            MultiBlockLattice2D<T,Descriptor>& lattice,
                            Box2D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computeVelocity (
                                  MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
                                  Box2D domain );

template<typename T, template<typename U> class Descriptor>
void computePiNeq(MultiBlockLattice2D<T,Descriptor>& lattice,
                             MultiNTensorField2D<T>& PiNeq, Box2D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputePiNeq(MultiBlockLattice2D<T,Descriptor>& lattice,
                                   MultiNTensorField2D<T>& PiNeq, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computePiNeq (
                            MultiBlockLattice2D<T,Descriptor>& lattice,
                            Box2D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computePiNeq (
                                  MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
                                  Box2D domain );


template<typename T, template<typename U> class Descriptor>
void computeShearStress(MultiBlockLattice2D<T,Descriptor>& lattice,
                  MultiNTensorField2D<T>& ShearStress, Box2D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputeShearStress(MultiBlockLattice2D<T,Descriptor>& lattice,
                        MultiNTensorField2D<T>& ShearStress, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computeShearStress (
    MultiBlockLattice2D<T,Descriptor>& lattice,
    Box2D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computeShearStress (
    MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
    Box2D domain );

template<typename T, template<typename U> class Descriptor>
void computeStrainRateFromStress(MultiBlockLattice2D<T,Descriptor>& lattice,
                                 MultiNTensorField2D<T>& S, Box2D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputeStrainRateFromStress(MultiBlockLattice2D<T,Descriptor>& lattice,
                                       MultiNTensorField2D<T>& S, MultiNTensorField2D<int>& mask, Box2D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computeStrainRateFromStress (
                           MultiBlockLattice2D<T,Descriptor>& lattice,
                           Box2D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computeStrainRateFromStress (
                                 MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
                                 Box2D domain );

template<typename T, template<typename U> class Descriptor>
void computePopulation(MultiBlockLattice2D<T,Descriptor>& lattice,
                       MultiNTensorField2D<T>& population,
                       Box2D domain, plint iPop);

template<typename T, template<typename U> class Descriptor>
void maskedComputePopulation(MultiBlockLattice2D<T,Descriptor>& lattice,
                             MultiNTensorField2D<T>& population, MultiNTensorField2D<int>& mask,
                             Box2D domain, plint iPop);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computePopulation (
                            MultiBlockLattice2D<T,Descriptor>& lattice,
                            Box2D domain, plint iPop );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computePopulation (
                                  MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
                                  Box2D domain, plint iPop );


template<typename T, template<typename U> class Descriptor>
void computePopulations(MultiBlockLattice2D<T,Descriptor>& lattice,
                        MultiNTensorField2D<T>& populations,
                        Box2D domain);

template<typename T, template<typename U> class Descriptor>
void maskedComputePopulations(MultiBlockLattice2D<T,Descriptor>& lattice,
                              MultiNTensorField2D<T>& populations, MultiNTensorField2D<int>& mask,
                              Box2D domain);

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computePopulations (
                            MultiBlockLattice2D<T,Descriptor>& lattice,
                            Box2D domain );

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computePopulations (
                                  MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
                                  Box2D domain );

}  // namespace plb

#endif  // SWIG_LATTICE_ANALYSIS_WRAPPER_2D_H
