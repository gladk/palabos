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
  * This program demonstrates how to cut surfaces represented as
  * sets of triangles by using cutting planes and cuboids.
  **/

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;

typedef double T;


int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./");

    plint n = 80; // Resolution (cylinder diameter).
    T radius = (T) n / 2.; // Cylinder radius.
    Array<T,3> center(0.0, 0.0, 0.0); // Cylinder center.
    T length = (T) (4 * n); // Cylinder length.
    plint nAxial = n / 2; // Two parameters needed for the creation of the triangularized cylinder surface.
    plint nCirc  = 3 * n / 2;

    // Create a cylinder surface as a set of triangles.
    TriangleSet<T> cylinder;
    cylinder = constructCylinder<T>(center, radius, radius, length, nAxial, nCirc);
    cylinder.writeAsciiSTL("cylinder.stl");

    // Define the cut plane by a point and a normal vector.
    Array<T,3> cutPlanePoint(length / 3.3, 0.0, 0.0);
    Array<T,3> cutPlaneNormal(-1.0, -0.1, -0.5);

    Plane<T> cutPlane;
    cutPlane.point = cutPlanePoint;
    cutPlane.normal = cutPlaneNormal;

    // Cut the cylinder in two parts.
    TriangleSet<T> rightCylinderPart;
    int rv = 0;
    rv = cylinder.cutWithPlane(cutPlane, rightCylinderPart);
    if (rv != 1) {
        pcout << "Problem with the surface cutting." << std::endl;
        exit(1);
    }
    rightCylinderPart.writeAsciiSTL("rightCylinderPart.stl");

    TriangleSet<T> leftCylinderPart;
    cutPlane.normal = -cutPlane.normal;
    rv = 0;
    rv = cylinder.cutWithPlane(cutPlane, leftCylinderPart);
    if (rv != 1) {
        pcout << "Problem with the surface cutting." << std::endl;
        exit(1);
    }
    leftCylinderPart.writeAsciiSTL("leftCylinderPart.stl");

    // Use a cuboid to cut the cylinder (correct usage).
    Array<T,3> lowerLeftCorner = cutPlanePoint - Array<T,3>(length, 5.0 * radius / 4.0, 5.0 * radius / 4.0);
    Array<T,3> upperRightCorner = cutPlanePoint + Array<T,3>(length / 4.0, 5.0 * radius / 4.0, 5.0 * radius / 4.0);

    Cuboid<T> cutCuboid;
    cutCuboid.lowerLeftCorner = lowerLeftCorner;
    cutCuboid.upperRightCorner = upperRightCorner;

    TriangleSet<T> newRightCylinderPart;
    cutPlane.normal = -cutPlane.normal;
    rv = 0;
    rv = cylinder.cutWithPlane(cutPlane, cutCuboid, newRightCylinderPart);
    if (rv != 1) {
        pcout << "Problem with the surface cutting." << std::endl;
        exit(1);
    }
    newRightCylinderPart.writeAsciiSTL("newRightCylinderPart.stl");

    // Use a cuboid to cut the cylinder (wrong usage).
    upperRightCorner = cutPlanePoint + Array<T,3>(0.0, radius, radius);

    cutCuboid.upperRightCorner = upperRightCorner;

    TriangleSet<T> newBrokenRightCylinderPart;
    rv = 0;
    rv = cylinder.cutWithPlane(cutPlane, cutCuboid, newBrokenRightCylinderPart);
    if (rv != 1) {
        pcout << "Problem with the surface cutting." << std::endl;
        exit(1);
    }
    newBrokenRightCylinderPart.writeAsciiSTL("newBrokenRightCylinderPart.stl");

    return 0;
}

