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

#ifndef POISEUILLE_H

#include "palabos2D.h"

/// Velocity on the parabolic Poiseuille profile
template<typename T>
T poiseuilleVelocity(int iY, plb::IncomprFlowParam<T> const& parameters);

/// Linearly decreasing pressure profile
template<typename T>
T poiseuillePressure(int iX, plb::IncomprFlowParam<T> const& parameters);

/// Convert pressure to density according to ideal gas law
template<typename T, template<typename U> class Descriptor>
T poiseuilleDensity(int iX, plb::IncomprFlowParam<T> const& parameters);

/// A functional, used to initialize the velocity for the boundary conditions
template<typename T>
class PoiseuilleVelocity;

/// A functional, used to initialize the density for the boundary conditions
template<typename T, template<typename U> class Descriptor>
class PoiseuilleDensity;

/// A functional, used to create an initial condition for the density and velocity
template<typename T, template<typename U> class Descriptor>
class PoiseuilleVelocityAndDensity;

template<typename T, template<typename U> class Descriptor>
void createPoiseuilleBoundaries( plb::MultiBlockLattice2D<T,Descriptor>& lattice,
                                 plb::IncomprFlowParam<T> const& parameters,
                                 plb::OnLatticeBoundaryCondition2D<T,Descriptor>& boundaryCondition );

template<typename T, template<typename U> class Descriptor>
void createPoiseuilleInitialValues( plb::MultiBlockLattice2D<T,Descriptor>& lattice,
                                    plb::IncomprFlowParam<T> const& parameters );

#endif  // POISEUILLE_H
