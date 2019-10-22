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

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./");
    global::IOpolicy().activateParallelIO(false);

    //std::cout.precision(16);
    //std::scientific(std::cout);

    std::string precisionStr;
    try {
        global::argv(1).read(precisionStr);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " [FLT | DBL | LDBL | INF] inputSTL.stl" << std::endl;
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

    std::string stlFileName;
    try {
        global::argv(2).read(stlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " [FLT | DBL | LDBL | INF] inputSTL.stl" << std::endl;
        exit(-1);
    }

    TriangleSet<T>* set = 0;
    try {
        set = new TriangleSet<T>(stlFileName, precision);
    }
    catch (PlbIOException& exception) {
        pcout << "ERROR, could not read STL file " << stlFileName
              << ": " << exception.what() << std::endl;
        exit(-1);
    }

    plint numTriangles = set->getTriangles().size();
    pcout << "The STL file contains " << numTriangles << " triangles." << std::endl;

    plint numZeroAreaTriangles = set->numZeroAreaTriangles();
    if (numZeroAreaTriangles) {
        pcout << "The TriangleSet has " << numZeroAreaTriangles << " zero-area triangles!" << std::endl;
    } else {
        pcout << "The TriangleSet does not have any zero-area triangles." << std::endl;
    }

    if (set->hasFloatingPointPrecisionDependence()) {
        pcout << "The triangle set has a dependence on the floating point precision used." << std::endl;
    } else {
        pcout << "The triangle set does not have a dependence on the floating point precision used." << std::endl;
    }

    pcout << "The minimum triangle edge length is: " << set->getMinEdgeLength() << std::endl;
    pcout << "The maximum triangle edge length is: " << set->getMaxEdgeLength() << std::endl;

    pcout << "The minimum triangle area is: " << set->getMinTriangleArea() << std::endl;
    pcout << "The maximum triangle area is: " << set->getMaxTriangleArea() << std::endl;

    Cuboid<T> bCuboid = set->getBoundingCuboid();
    Array<T,3> llc = bCuboid.lowerLeftCorner;
    Array<T,3> urc = bCuboid.upperRightCorner;
    pcout << "The bounding box of the triangle set is: " << std::endl;
    pcout << "  x: [" << llc[0] << ", " << urc[0] << "]" << std::endl;
    pcout << "  y: [" << llc[1] << ", " << urc[1] << "]" << std::endl;
    pcout << "  z: [" << llc[2] << ", " << urc[2] << "]" << std::endl;
    Array<T,3> c = (T) 0.5 * (llc + urc);
    pcout << "The center of the bounding box is: [" << c[0] << ", " << c[1] << ", " << c[2] << "]" << std::endl;

    c = set->getCentroid();
    pcout << "The centroid of the triangle set is: [" << c[0] << ", " << c[1] << ", " << c[2] << "]" << std::endl;

    delete set;

    return 0;
}

