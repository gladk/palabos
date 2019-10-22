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
 * Division of a plane for interpolation in 3D grid refinement -- header file.
 */
#ifndef DOMAIN_DIVISION_3D_H
#define DOMAIN_DIVISION_3D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/geometry3D.h"

namespace plb {

/// Fragmentation of the surface of a block into bulk, edges, and corners.
class SurfaceDivision3D {
public:
    SurfaceDivision3D(Box3D const& domain_, plint boundaryWidth, plint direction_);
    Box3D bulk() const;
    Box3D edge0N() const;
    Box3D edge0P() const;
    Box3D edge1N() const;
    Box3D edge1P() const;
    Box3D cornerNN() const;
    Box3D cornerPN() const;
    Box3D cornerNP() const;
    Box3D cornerPP() const;
private:
    Box3D domain;
    plint bw;
    plint direction;
};

}  // namespace plb

#endif //  DOMAIN_DIVISION_3D_H
