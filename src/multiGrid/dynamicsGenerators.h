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

/* Main author: Daniel Lagrava
 **/

/** \file
 * Generator of Dynamics vectors for multigrid 
 */
#ifndef DYNAMICS_GENERATORS_H
#define DYNAMICS_GENERATORS_H

#include "basicDynamics/isoThermalDynamics.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
std::vector<Dynamics<T,Descriptor>*> generateBGKdynamics(RefinementParameters<T>* refinementParameters){
    
    plint numLevel = refinementParameters->getNumLevels();
    
    std::vector<Dynamics<T,Descriptor>*> result(numLevel);
    for (plint iLevel=0; iLevel<numLevel; ++iLevel) {
        T omega = refinementParameters->getParameters(iLevel).getOmega();
        result[iLevel] = new BGKdynamics<T,Descriptor>(omega);
    }
    
    return result;
}


template <typename T, template <typename U> class Descriptor>
std::vector<Dynamics<T,Descriptor>*> generateBounceBackDynamics(RefinementParameters<T>* refinementParameters){
    
    plint numLevel = refinementParameters->getNumLevels();
    
    std::vector<Dynamics<T,Descriptor>*> result(numLevel);
    for (plint iLevel=0; iLevel<numLevel; ++iLevel) {
        result[iLevel] = new BounceBack<T,Descriptor>(0.0);        
    }
    
    return result;
}


//TODO add generators for all dynamics that will be used

} // namespace plb

#endif // DYNAMICS_GENERATORS_H
