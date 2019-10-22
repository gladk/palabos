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
#ifndef SWIG_LATTICE_INITIALIZER_WRAPPER_2D_H
#define SWIG_LATTICE_INITIALIZER_WRAPPER_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "core/dynamics.h"
#include "plbWrapper/lattice/latticeInitializerFunctional2D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
void pypalDefineDynamics( MultiBlockLattice2D<T,Descriptor>& lattice,
                          Box2D domain, Dynamics<T,Descriptor>* dynamics );

template<typename T, template<typename U> class Descriptor>
void maskedDefineDynamics( MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
                           Box2D domain, Dynamics<T,Descriptor>* dynamics );

template<typename T, template<class U> class Descriptor>
void setBoundaryVelocity( MultiBlockLattice2D<T,Descriptor>& lattice,
                          T* velocity, int numDimIs2, Box2D domain );

template<typename T, template<class U> class Descriptor>
void setBoundaryVelocity( MultiBlockLattice2D<T,Descriptor>& lattice,
                          MultiNTensorField2D<T>& velocity,
                          Box2D domain );

template<typename T, template<class U> class Descriptor>
void maskedSetBoundaryVelocity( MultiBlockLattice2D<T,Descriptor>& lattice,
                                MultiNTensorField2D<T>& velocity,
                                MultiNTensorField2D<int>& mask,
                                Box2D domain );

template<typename T, template<class U> class Descriptor>
void initializeAtEquilibrium( MultiBlockLattice2D<T,Descriptor>& lattice,
                              T density, T* velocity, int numDimIs2, Box2D domain );

template<typename T, template<class U> class Descriptor>
void initializeAtEquilibrium( MultiBlockLattice2D<T,Descriptor>& lattice,
                              MultiNTensorField2D<T>& density,
                              MultiNTensorField2D<T>& velocity,
                              Box2D domain );

template<typename T, template<class U> class Descriptor>
void maskedInitializeAtEquilibrium( MultiBlockLattice2D<T,Descriptor>& lattice,
                                    MultiNTensorField2D<T>& density,
                                    MultiNTensorField2D<T>& velocity,
                                    MultiNTensorField2D<int>& mask,
                                    Box2D domain );

template<typename T, template<class U> class Descriptor>
void setExternalVector( MultiBlockLattice2D<T,Descriptor>& lattice,
                        int vectorStartsAt, T* externalVector, int numDimIs2, Box2D domain );

template<typename T, template<class U> class Descriptor>
void setPopulations( MultiBlockLattice2D<T,Descriptor>& lattice,
                     T* populations, int numDimIsQ, Box2D domain );

template<typename T, template<class U> class Descriptor>
void setPopulations( MultiBlockLattice2D<T,Descriptor>& lattice,
                     MultiNTensorField2D<T>& populations,
                     Box2D domain );

template<typename T, template<class U> class Descriptor>
void maskedSetPopulations( MultiBlockLattice2D<T,Descriptor>& lattice,
                           T* populations, int numDimIsQ,
                           MultiNTensorField2D<int>& mask, Box2D domain );

template<typename T, template<class U> class Descriptor>
void maskedSetPopulations( MultiBlockLattice2D<T,Descriptor>& lattice,
                           MultiNTensorField2D<T>& populations,
                           MultiNTensorField2D<int>& mask,
                           Box2D domain );

}  // namespace plb

#endif  // SWIG_LATTICE_INITIALIZER_WRAPPER_2D_H
