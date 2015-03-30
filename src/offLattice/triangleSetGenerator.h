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

/* Main author: Dimitrios Kontaxakis */

#ifndef TRIANGLE_SET_GENERATOR_H
#define TRIANGLE_SET_GENERATOR_H

#include "core/globalDefs.h"
#include "core/array.h"
#include "offLattice/triangleSet.h"

#include <vector>

namespace plb {

/// Create and return a sphere as a set of triangles. The center and radius of the sphere
///   must be given as arguments. A minimum number of triangles for the sphere triangulation
///   must be provided as well. This number is suggestive for the resolution. The
///   actual number of triangles can be greater than the one provided.
template<typename T>
TriangleSet<T> constructSphere(Array<T,3> const& center, T radius, plint minNumOfTriangles);

/// Create and return a circular tapered cylinder. The axis of the cylinder is parallel
///   to the x-axis. The center of the inlet disk must be given as argument, as well as
///   the inlet and outer radii and the length of the object. ``nAxial'' is the
///   number of points in the axial direction before triangulation, and ``nCirc''
///   is the number of points in the circumferential direction before triangulation.
///   The final number of points will be 2*nAxial-1 and 2*nCirc in the axial
///   and circumferential directions, respectively.
template<typename T>
TriangleSet<T> constructCylinder(Array<T,3> const& inletCenter, T inletRadius, T outletRadius,
                                 T length, plint nAxial, plint nCirc);

template<typename T>
TriangleSet<T> constructCylinder( Array<T,3> const& inletCenter, T inletRadius, T outletRadius,
                                  T length, plint nAxial, plint nCirc,
                                  std::vector<Array<T,3> >& inletPoints );

template<typename T>
TriangleSet<T> constructCylinder( Array<T,3> const& inletCenter, Array<T,3> const& axis,
                                  T inletRadius, T outletRadius,
                                  T length, plint nAxial, plint nCirc,
                                  std::vector<Array<T,3> >& inletPoints );

template<typename T>
TriangleSet<T> constructCuboid (
        Array<T,3> const& lowerCorner, Array<T,3> const& upperCorner,
        Array<plint,3> const& nSegments );

template<typename T>
TriangleSet<T> patchTubes(TriangleSet<T> const& geometryWithOpenings, plint sortDirection, std::vector<T> patchLengths);

/// Create and return a rectangle. The rectangle is on the x-y plane, and its lower left
///   corner is at the origin of the axes. Its sides have length "lx" and "ly", while
///   the number of points for the triangulation are "nx" and "ny" on the x and y axis,
///   respectively. This means that the total number of triangles is 2*(nx-1)*(ny-1).
template<typename T>
TriangleSet<T> constructRectangle(T lx, T ly, plint nx, plint ny);

/// Create and return a generically placed rectangle. The center of the rectangle is placed
///   at "center", and its normal is "normal". Its sides have length "lx" and "ly", while
///   the number of points for the triangulation are "nx" and "ny" on the x and y axis,
///   respectively. This means that the total number of triangles is 2*(nx-1)*(ny-1).
template<typename T>
TriangleSet<T> constructGenericRectangle(T lx, T ly, plint nx, plint ny, Array<T,3> const& center, Array<T,3> const& normal);

/// Create a strip of triangles. There are given two sets of points "from" and "to". These
///   must contain the same number of points. This function creates a set of triangles
///   by connecting points from "from" with points from "to".
template<typename T>
TriangleSet<T> constructStrip(std::vector<Array<T,3> > from, std::vector<Array<T,3> > to);

} // namespace plb

#endif  // TRIANGLE_SET_GENERATOR_H

