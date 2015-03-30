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

#ifndef CREATE_BUBBLES_3D_H
#define CREATE_BUBBLES_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include <limits>


namespace plb {

template<typename T, template<typename U> class Descriptor>
void punchSphere( std::vector<MultiBlock3D*> const& twoPhaseArgs, Array<T,3> const& center, T radius,
                  T rhoEmpty, T rho0, Dynamics<T,Descriptor>& dynamics, Box3D domain );

template<typename T, template<typename U> class Descriptor>
void analyticalPunchSphere( std::vector<MultiBlock3D*> const& twoPhaseArgs, Array<T,3> const& center, T radius,
                            T rhoEmpty, T rho0, plint subDivision, Dynamics<T,Descriptor>& dynamics, Box3D domain );

template<typename T, template<typename U> class Descriptor>
T computeAverageSphereDensity( std::vector<MultiBlock3D*> const& twoPhaseArgs,
                               Array<T,3> const& center, T radius, Box3D domain );

}  // namespace plb

#endif  // CREATE_BUBBLES_3D_H

