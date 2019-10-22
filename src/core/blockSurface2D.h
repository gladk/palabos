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
 * Handling the surface area of a block -- header file.
 */
#ifndef BLOCK_SURFACE_2D_H
#define BLOCK_SURFACE_2D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/geometry2D.h"

namespace plb {

/// Fragmentation of the surface of a block into bulk, edges, and corners.
class BlockSurface2D {
public:
    BlockSurface2D(Box2D const& domain_, plint boundaryWidth);
    Box2D bulk() const;
    Box2D edge0N() const;
    Box2D edge0P() const;
    Box2D edge1N() const;
    Box2D edge1P() const;
    Box2D cornerNN() const;
    Box2D cornerPN() const;
    Box2D cornerNP() const;
    Box2D cornerPP() const;
private:
    Box2D domain;
    plint bw;
};

}  // namespace plb

#endif //  BLOCK_SURFACE_2D_H
