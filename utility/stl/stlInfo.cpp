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

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;

typedef double T;

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./");
    global::IOpolicy().activateParallelIO(false);

    std::cout.precision(16);
    std::scientific(std::cout);

    std::string stlFileName;
    try {
        global::argv(1).read(stlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " inputSTL.stl" << std::endl;
        exit(-1);
    }

    TriangleSet<T>* set = 0;
    try {
        set = new TriangleSet<T>(stlFileName, DBL);
    }
    catch (PlbIOException& exception) {
        pcout << "ERROR, could not read STL file " << stlFileName
              << ": " << exception.what() << std::endl;
        exit(-1);
    }

    plint numTriangles = set->getTriangles().size();
    pcout << "The STL file contains " << numTriangles << " triangles." << std::endl;

    if (set->hasDegenerateTriangles()) {
        pcout << "The triangle set contains degenerate triangles." << std::endl;
    } else {
        pcout << "The triangle set does not contain degenerate triangles." << std::endl;
    }

    pcout << "The minimum triangle edge length is: " << set->getMinEdgeLength() << std::endl;
    pcout << "The maximum triangle edge length is: " << set->getMaxEdgeLength() << std::endl;

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

