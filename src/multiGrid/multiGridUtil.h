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
 * Set of functions commonly used in LB computations
 *  -- header file
 */
#ifndef MULTI_GRID_UTIL_H
#define MULTI_GRID_UTIL_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/util.h"

namespace plb {

/// dxScale is positive if the system is coarser than the reference.
inline double scaleToReference(int dxScale, int dimDx, int dtScale, int dimDt) {
    return util::twoToThePower(dimDx*dxScale+dimDt*dtScale);
}

/// dxScale is positive if the system is coarser than the reference.
inline double scaleFromReference(int dxScale, int dimDx, int dtScale, int dimDt) {
    return util::twoToThePower(-dimDx*dxScale-dimDt*dtScale);
}

}  // namespace plb

#endif  // MULTI_GRID_UTIL_H
