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

/** \file Constants and definitions for boundary conditions -- header file.  */

#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

namespace plb {
   
    namespace boundary {
        /** dirichlet:     impose a velocity (* obsolete ...or a density *).
         *  neumann:       zero-gradient for all velocity components or for the density.
         *  freeslip:      zero-gradient for tangential velocity components, and zero for normal ones.
         *  density:       fixed density, Neumann for tangential velocity components.
         *  (*obsolete*) outflow:       zero-gradient for all velocity components.
         *  (*obsolete*) normalOutflow: zero-gradient for normal velocity components, and zero for tangential ones.
         **/
        typedef enum {dirichlet, neumann, freeslip, density, outflow, normalOutflow} BcType;
    }

}  // namespace plb

#endif
