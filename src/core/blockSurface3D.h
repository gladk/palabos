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
#ifndef BLOCK_SURFACE_3D_H
#define BLOCK_SURFACE_3D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/geometry3D.h"

namespace plb {

/// Fragmentation of the surface of a block into bulk, surfaces, edges, and corners.
class BlockSurface3D {
public:
    BlockSurface3D(Box3D const& domain_, plint boundaryWidth);
    Box3D bulk() const;
    Box3D surface0N() const;
    Box3D surface0P() const;
    Box3D surface1N() const;
    Box3D surface1P() const;
    Box3D surface2N() const;
    Box3D surface2P() const;
    Box3D edge0NN() const;
    Box3D edge0NP() const;
    Box3D edge0PN() const;
    Box3D edge0PP() const;
    Box3D edge1NN() const;
    Box3D edge1NP() const;
    Box3D edge1PN() const;
    Box3D edge1PP() const;
    Box3D edge2NN() const;
    Box3D edge2NP() const;
    Box3D edge2PN() const;
    Box3D edge2PP() const;
    Box3D cornerNNN() const;
    Box3D cornerNNP() const;
    Box3D cornerNPN() const;
    Box3D cornerNPP() const;
    Box3D cornerPNN() const;
    Box3D cornerPNP() const;
    Box3D cornerPPN() const;
    Box3D cornerPPP() const;
private:
    Box3D domain;
    plint bw;
};

}  // namespace plb

#endif //  BLOCK_SURFACE_3D_H
