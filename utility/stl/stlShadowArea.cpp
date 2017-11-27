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

#include "palabos3D.h"
#include "palabos3D.hh"

#define AWAY_FROM_SURFACE   0
#define CLOSE_TO_SURFACE    1

using namespace plb;

typedef double T;

typedef TriangleSet<T>::Triangle Triangle;

T toPhys(T lbVal, plint direction, T dx, Array<T,3> const& location)
{
    PLB_ASSERT(direction >= 0 && direction <= 2);
    return (lbVal * dx + location[direction]);
}

Array<T,3> toPhys(Array<T,3> const& lbVal, T dx, Array<T,3> const& location)
{
    return (lbVal * dx + location);
}

T toLB(T physVal, plint direction, T dx, Array<T,3> const& location)
{
    PLB_ASSERT(direction >= 0 && direction <= 2);
    return (physVal - location[direction]) / dx;
}

Array<T,3> toLB(Array<T,3> const& physVal, T dx, Array<T,3> const& location)
{
    return (physVal - location) / dx;
}

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./");
    global::IOpolicy().activateParallelIO(false);

    if (argc < 5) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " [FLT | DBL | LDBL | INF] axis N stl1.stl stl2.stl ..." << std::endl;
        pcout << "The parameter axis is 0, 1, or 2, and denotes the axis normal to the plane for the projection." << std::endl;
        pcout << "The parameter N is the resolution of the minimum bounding box length." << std::endl;
        pcout << "N should be very high (~1000)." << std::endl;
        exit(-1);
    }

    std::string precisionStr;
    try {
        global::argv(1).read(precisionStr);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " [FLT | DBL | LDBL | INF] axis N stl1.stl stl2.stl ..." << std::endl;
        pcout << "The parameter axis is 0, 1, or 2, and denotes the axis normal to the plane for the projection." << std::endl;
        pcout << "The parameter N is the resolution of the minimum bounding box length." << std::endl;
        pcout << "N should be very high (~1000)." << std::endl;
        exit(-1);
    }

    Precision precision;
    if (precisionStr == "FLT") {
        precision = FLT;
    } else if (precisionStr == "DBL") {
        precision = DBL;
    } else if (precisionStr == "LDBL") {
        precision = LDBL;
    } else if (precisionStr == "INF") {
        precision = INF;
    } else {
        pcout << "Wrong precision command-line argument." << std::endl;
        exit(-1);
    }

    plint axis = -1;
    global::argv(2).read(axis);
    PLB_ASSERT(axis == 0 || axis == 1 || axis == 2);

    plint N = -1;
    global::argv(3).read(N);
    PLB_ASSERT(N > 2);

    TriangleSet<T>* set = new TriangleSet<T>(precision);
    for (int i = 4; i < argc; i++) {
        pcout << "Reading " << argv[i] << std::endl;
        TriangleSet<T>* mesh = new TriangleSet<T>(argv[i], precision);
        set->append(*mesh);
        delete mesh;
    }

    pcout << "Computing the projected geometry." << std::endl;
    std::vector<Triangle> const& triangles = set->getTriangles();
    std::vector<Triangle> projectedTriangles(triangles.size());
    for (plint iTriangle = 0; iTriangle < (plint) triangles.size(); iTriangle++) {
        Triangle t(triangles[iTriangle]);
        t[0][axis] = 0.0;
        t[1][axis] = 0.0;
        t[2][axis] = 0.0;
        projectedTriangles[iTriangle] = t;
    }

    delete set; set = 0;

    TriangleSet<T>* projectedSet = new TriangleSet<T>(projectedTriangles, precision);
    pcout << "Saving projectedGeometry.stl" << std::endl;
    projectedSet->writeAsciiSTL("projectedGeometry.stl");

    Cuboid<T> cuboid = projectedSet->getBoundingCuboid();
    T lx = cuboid.x1() - cuboid.x0();
    T ly = cuboid.y1() - cuboid.y0();
    T lz = cuboid.z1() - cuboid.z0();
    PLB_ASSERT((axis == 0 && util::isZero(lx)) || (axis == 1 && util::isZero(ly)) || (axis == 2 && util::isZero(lz)));
    T delta = -1.0;
    if (axis == 0) {
        delta = std::min(ly, lz);
    } else if (axis == 1) {
        delta = std::min(lx, lz);
    } else {
        delta = std::min(lx, ly);
    }

    T scalingFactor = (T) N / delta;
    Array<T,3> physicalLocation = cuboid.lowerLeftCorner;
    T dx = (T) 1 / scalingFactor;
    plint symmetricLayer = 2;
    Array<T,3> luOffset(symmetricLayer, symmetricLayer, symmetricLayer);
    luOffset[axis] = 0.0;
    projectedSet->translate(-physicalLocation);
    projectedSet->scale(scalingFactor);
    projectedSet->translate(luOffset);
    physicalLocation -= luOffset * dx;

    Cuboid<T> cuboid_LB;
    cuboid_LB.lowerLeftCorner = toLB(cuboid.lowerLeftCorner, dx, physicalLocation);
    cuboid_LB.upperRightCorner = toLB(cuboid.upperRightCorner, dx, physicalLocation);

    T dX = cuboid_LB.x1() - cuboid_LB.x0();
    T dY = cuboid_LB.y1() - cuboid_LB.y0();
    T dZ = cuboid_LB.z1() - cuboid_LB.z0();

    plint nx = (plint) dX + 1 + (axis != 0 ? 2 * symmetricLayer : 0);
    plint ny = (plint) dY + 1 + (axis != 1 ? 2 * symmetricLayer : 0);
    plint nz = (plint) dZ + 1 + (axis != 2 ? 2 * symmetricLayer : 0);

    pcout << "Resolution for the area computation: " << dx << std::endl;

    pcout << "Refining the geometry." << std::endl;
    plint maxNumRef = 100;
    T dx_LB = 1.0;
    bool succeeded = projectedSet->refineRecursively(dx_LB, maxNumRef);
    if (!succeeded) {
        pcerr << "Error: The target maximum triangle edge length " << dx
              << "       in lattice units was not reached after " << maxNumRef << " refinement iterations." << std::endl;
        exit(1);
    }

    pcout << "Computing the connected set." << std::endl;
    ConnectedTriangleSet<T>* connectedTriangleSet = new ConnectedTriangleSet<T>(*projectedSet);
    plint numVertices = connectedTriangleSet->getNumVertices();
    std::vector<Array<T,3> > vertices(numVertices);
    for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
        vertices[iVertex] = connectedTriangleSet->getVertex(iVertex);
    }
    delete connectedTriangleSet; connectedTriangleSet = 0;

    MultiScalarField3D<int>* nodeTags = new MultiScalarField3D<int>(nx, ny, nz, AWAY_FROM_SURFACE);
    nodeTags->toggleInternalStatistics(false);

    MultiContainerBlock3D* container = new MultiContainerBlock3D(*nodeTags);
    container->toggleInternalStatistics(false);

    pcout << "Creating the integer tags." << std::endl;
    plint halfWidth = 1;
    instantiateSurfaceBlockData<T>(halfWidth, vertices, *container);
    surfaceOnLattice<T,int>(CLOSE_TO_SURFACE, halfWidth, *nodeTags, *container);
    delete container; container = 0;

    pcout << "Saving projectedGeometry.vti" << std::endl;
    {
        VtkImageOutput3D<T> vtkOut("projectedGeometry", dx, physicalLocation);
        vtkOut.writeData<float>(*copyConvert<int,T>(*nodeTags, nodeTags->getBoundingBox()), "flag", 1.0);
    }

    pcout << std::endl;
    T sum = computeSum(*nodeTags);
    pcout << "Projected area = " << sum * dx * dx << std::endl;

    delete nodeTags; nodeTags = 0;

    return 0;
}
