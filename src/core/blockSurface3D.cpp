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
 * Handling the surface area of a block -- implementation.
 */

#include "core/blockSurface3D.h"

namespace plb {

BlockSurface3D::BlockSurface3D(Box3D const& domain_, plint boundaryWidth)
    : domain(domain_),
      bw(boundaryWidth)
{ }

Box3D BlockSurface3D::bulk() const {
    return Box3D( domain.x0+bw,   domain.x1-bw,
                  domain.y0+bw,   domain.y1-bw,
                  domain.z0+bw,   domain.z1-bw    );
}

Box3D BlockSurface3D::surface0N() const {
    return Box3D( domain.x0,      domain.x0+bw-1,
                  domain.y0+bw,   domain.y1-bw,
                  domain.z0+bw,   domain.z1-bw    );
}

Box3D BlockSurface3D::surface0P() const {
    return Box3D( domain.x1-bw+1, domain.x1,
                  domain.y0+bw,   domain.y1-bw,
                  domain.z0+bw,   domain.z1-bw    );
}

Box3D BlockSurface3D::surface1N() const {
    return Box3D( domain.x0+bw,   domain.x1-bw,
                  domain.y0,      domain.y0+bw-1,
                  domain.z0+bw,   domain.z1-bw    );
}

Box3D BlockSurface3D::surface1P() const {
    return Box3D( domain.x0+bw,   domain.x1-bw,
                  domain.y1-bw+1, domain.y1,
                  domain.z0+bw,   domain.z1-bw    );
}

Box3D BlockSurface3D::surface2N() const {
    return Box3D( domain.x0+bw,   domain.x1-bw,
                  domain.y0+bw,   domain.y1-bw,
                  domain.z0,      domain.z0+bw-1  );
}

Box3D BlockSurface3D::surface2P() const {
    return Box3D( domain.x0+bw,   domain.x1-bw,
                  domain.y0+bw,   domain.y1-bw,
                  domain.z1-bw+1, domain.z1       );
}

Box3D BlockSurface3D::edge0NN() const {
    return Box3D( domain.x0+bw,   domain.x1-bw,       // Plane x;  y: negative, z: negative
                  domain.y0,      domain.y0+bw-1,
                  domain.z0,      domain.z0+bw-1  );
}

Box3D BlockSurface3D::edge0NP() const {
    return Box3D( domain.x0+bw,   domain.x1-bw,       // Plane x;  y: negative, z: positive
                  domain.y0,      domain.y0+bw-1,
                  domain.z1-bw+1, domain.z1       );
}

Box3D BlockSurface3D::edge0PN() const {
    return Box3D( domain.x0+bw,   domain.x1-bw,       // Plane x;  y: positive, z: negative
                  domain.y1-bw+1, domain.y1,
                  domain.z0,      domain.z0+bw-1  );
}

Box3D BlockSurface3D::edge0PP() const {
    return Box3D( domain.x0+bw,   domain.x1-bw,       // Plane x;  y: positive, z: positive
                  domain.y1-bw+1, domain.y1,
                  domain.z1-bw+1, domain.z1       );
}

Box3D BlockSurface3D::edge1NN() const {
    return Box3D( domain.x0,      domain.x0+bw-1,     // Plane y;  z: negative, x: negative
                  domain.y0+bw,   domain.y1-bw,
                  domain.z0,      domain.z0+bw-1  );
}

Box3D BlockSurface3D::edge1NP() const {
    return Box3D( domain.x1-bw+1, domain.x1,          // Plane y;  z: negative, x: positive
                  domain.y0+bw,   domain.y1-bw,
                  domain.z0,      domain.z0+bw-1  );
}

Box3D BlockSurface3D::edge1PN() const {
    return Box3D( domain.x0,      domain.x0+bw-1,     // Plane y;  z: positive, x: negative
                  domain.y0+bw,   domain.y1-bw,
                  domain.z1-bw+1, domain.z1       );
}

Box3D BlockSurface3D::edge1PP() const {
    return Box3D( domain.x1-bw+1, domain.x1,          // Plane y;  z: positive, x: positive
                  domain.y0+bw,   domain.y1-bw,
                  domain.z1-bw+1, domain.z1       );
}

Box3D BlockSurface3D::edge2NN() const {
    return Box3D( domain.x0,      domain.x0+bw-1,     // Plane z;  x: negative, y: negative
                  domain.y0,      domain.y0+bw-1,
                  domain.z0+bw,   domain.z1-bw    );
}

Box3D BlockSurface3D::edge2NP() const {
    return Box3D( domain.x0,      domain.x0+bw-1,     // Plane z;  x: negative, y: positive
                  domain.y1-bw+1, domain.y1,
                  domain.z0+bw,   domain.z1-bw    );
}

Box3D BlockSurface3D::edge2PN() const {
    return Box3D( domain.x1-bw+1, domain.x1,          // Plane z;  x: positive, y: negative
                  domain.y0,      domain.y0+bw-1,
                  domain.z0+bw,   domain.z1-bw    );
}

Box3D BlockSurface3D::edge2PP() const {
    return Box3D( domain.x1-bw+1, domain.x1,          // Plane z;  x: positive, y: positive
                  domain.y1-bw+1, domain.y1,
                  domain.z0+bw,   domain.z1-bw    );
}

Box3D BlockSurface3D::cornerNNN() const {
    return Box3D( domain.x0,      domain.x0+bw-1,
                  domain.y0,      domain.y0+bw-1,
                  domain.z0,      domain.z0+bw-1  );
}

Box3D BlockSurface3D::cornerNNP() const {
    return Box3D( domain.x0,      domain.x0+bw-1,
                  domain.y0,      domain.y0+bw-1,
                  domain.z1-bw+1, domain.z1       );
}

Box3D BlockSurface3D::cornerNPN() const {
    return Box3D( domain.x0,      domain.x0+bw-1,
                  domain.y1-bw+1, domain.y1,
                  domain.z0,      domain.z0+bw-1  );
}

Box3D BlockSurface3D::cornerNPP() const {
    return Box3D( domain.x0,      domain.x0+bw-1,
                  domain.y1-bw+1, domain.y1,
                  domain.z1-bw+1, domain.z1       );
}

Box3D BlockSurface3D::cornerPNN() const {
    return Box3D( domain.x1-bw+1, domain.x1,
                  domain.y0,      domain.y0+bw-1,
                  domain.z0,      domain.z0+bw-1  );
}

Box3D BlockSurface3D::cornerPNP() const {
    return Box3D( domain.x1-bw+1, domain.x1,
                  domain.y0,      domain.y0+bw-1,
                  domain.z1-bw+1, domain.z1       );
}

Box3D BlockSurface3D::cornerPPN() const {
    return Box3D( domain.x1-bw+1, domain.x1,
                  domain.y1-bw+1, domain.y1,
                  domain.z0,      domain.z0+bw-1  );
}

Box3D BlockSurface3D::cornerPPP() const {
    return Box3D( domain.x1-bw+1, domain.x1,
                  domain.y1-bw+1, domain.y1,
                  domain.z1-bw+1, domain.z1       );
}

}  // namespace plb
