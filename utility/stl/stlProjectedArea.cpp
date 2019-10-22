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

#include <algorithm>
#include <string>

using namespace plb;

typedef double T;

typedef TriangleSet<T>::Triangle Triangle;

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./");
    global::IOpolicy().activateParallelIO(false);

    if (argc < 6) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " [FLT | DBL | LDBL | INF] nx ny nz stl1.stl stl2.stl ..." << std::endl;
        pcout << "The nx, ny, nz parameters are the components of the normal vector of the plane. " << std::endl;
        exit(-1);
    }

    std::string precisionStr;
    try {
        global::argv(1).read(precisionStr);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " [FLT | DBL | LDBL | INF] nx ny nz stl1.stl stl2.stl ..." << std::endl;
        pcout << "The nx, ny, nz parameters are the components of the normal vector of the plane. " << std::endl;
        exit(-1);
    }

    std::transform(precisionStr.begin(), precisionStr.end(), precisionStr.begin(), ::toupper);

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

    Array<T,3> planeNormal;
    global::argv(2).read(planeNormal[0]);
    global::argv(3).read(planeNormal[1]);
    global::argv(4).read(planeNormal[2]);
    T normNormal = norm(planeNormal);
    PLB_ASSERT(!util::isZero(normNormal));
    planeNormal /= normNormal;

    TriangleSet<T>* set = new TriangleSet<T>(precision);
    for (int i = 5; i < argc; i++) {
        TriangleSet<T>* mesh = new TriangleSet<T>(argv[i], precision);
        set->append(*mesh);
        delete mesh;
        pcout << "Read: " << argv[i] << std::endl;
    }

    plint numZeroAreaTriangles = set->numZeroAreaTriangles();
    if (numZeroAreaTriangles) {
        pcout << std::endl;
        pcout << "WARNING: the TriangleSet has " << numZeroAreaTriangles << " zero-area triangles!" << std::endl;
        pcout << std::endl;
    }

    if (set->hasFloatingPointPrecisionDependence()) {
        pcout << "WARNING: The TriangleSet has a dependence on the floating point precision used." << std::endl;
    }

    ConnectedTriangleSet<T>* connectedSet = 0;
    try {
        connectedSet = new ConnectedTriangleSet<T>(*set);
    }
    catch (PlbIOException& exception) {
        pcout << "ERROR: could not create a ConnectedTriangleSet object: "
              << exception.what() << std::endl;
        exit(-1);
    }
    delete set; set = 0;
    pcout << "The connected mesh has: " << connectedSet->getNumVertices() << " vertices, and "
          << connectedSet->getNumTriangles() << " triangles." << std::endl;

    T projectedArea = 0.0;
    for (plint iTriangle = 0; iTriangle < connectedSet->getNumTriangles(); iTriangle++) {
        T area;
        Array<T,3> unitNormal;
        connectedSet->computeTriangleAreaAndUnitNormal(iTriangle, area, unitNormal);
        projectedArea += dot(planeNormal, unitNormal) * area;
    }
    pcout << "The projected area is: " << projectedArea << std::endl;

    delete connectedSet;

    return 0;
}

