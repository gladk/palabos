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
 * Helper classes for parallel 3D multiblock lattice -- header file.
 */

#ifndef COMMUNICATION_PACKAGE_3D_H
#define COMMUNICATION_PACKAGE_3D_H

#include "core/geometry3D.h"
#include <vector>

namespace plb {

struct CommunicationInfo3D {
    plint fromBlockId;
    int fromProcessId;
    Box3D fromDomain;
    plint toBlockId;
    int toProcessId;
    Box3D toDomain;
    Dot3D absoluteOffset;
};

typedef std::vector<CommunicationInfo3D> CommunicationPackage3D;

}  // namespace plb

#endif  // COMMUNICATION_PACKAGE_3D_H
