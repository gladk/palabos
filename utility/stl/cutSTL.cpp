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

using namespace plb;

typedef double T;

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./");
    global::IOpolicy().activateParallelIO(false);

    if (argc != 9) {
        pcout << "Wrong parameters; the syntax is:" << std::endl
              << (std::string)global::argv(0)
              << " [FLT | DBL | LDBL | INF] inputSTL.stl cut-plane-point-x cut-plane-point-y cut-plane-point-z" << std::endl
              << "cut-plane-normal-x cut-plane-normal-y cut-plane-normal-z" << std::endl;
        exit(-1);
    }

    std::string precisionStr;
    std::string inStlFileName;
    Plane<T> cutPlane;
    try {
        global::argv(1).read(precisionStr);
        global::argv(2).read(inStlFileName);
        global::argv(3).read(cutPlane.point[0]);
        global::argv(4).read(cutPlane.point[1]);
        global::argv(5).read(cutPlane.point[2]);
        global::argv(6).read(cutPlane.normal[0]);
        global::argv(7).read(cutPlane.normal[1]);
        global::argv(8).read(cutPlane.normal[2]);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters." << std::endl;
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

    TriangleSet<T>* set = 0;
    try {
        set = new TriangleSet<T>(inStlFileName, precision);
    }
    catch (PlbIOException& exception) {
        pcout << "ERROR, could not read STL file " << inStlFileName
              << ": " << exception.what() << std::endl;
        exit(-1);
    }

    // Cut the geometry in two parts.
    TriangleSet<T> normalNegativePart;
    int rv = 0;
    rv = set->cutWithPlane(cutPlane, normalNegativePart);
    if (rv != 1) {
        pcout << "Problem with the surface cutting." << std::endl;
        exit(1);
    }
    normalNegativePart.writeAsciiSTL("normalNegativePart.stl");

    TriangleSet<T> normalPositivePart;
    cutPlane.normal = -cutPlane.normal;
    rv = 0;
    rv = set->cutWithPlane(cutPlane, normalPositivePart);
    if (rv != 1) {
        pcout << "Problem with the surface cutting." << std::endl;
        exit(1);
    }
    normalPositivePart.writeAsciiSTL("normalPositivePart.stl");

    delete set;

    return 0;
}

