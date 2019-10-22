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
using namespace std;

typedef double T;

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./");
    global::IOpolicy().activateParallelIO(false);

    string stlFileName, outFileName;
    try {
        global::argv(1).read(stlFileName);
        global::argv(2).read(outFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " inputSTL.stl outputSTL.stl [threshold [maxIter]]" << std::endl;
        exit(-1);
    }

    bool readThreshold = false;
    T threshold = 0.0;
    try {
        global::argv(3).read(threshold);
        readThreshold = true;
    }
    catch (PlbIOException& exception) {
        readThreshold = false;
    }

    bool readMaxIter = false;
    plint maxIter = 0;
    try {
        global::argv(4).read(maxIter);
        readMaxIter = true;
    }
    catch (PlbIOException& exception) {
        readMaxIter = false;
    }

    TriangleSet<T>* triangleSet = 0;
    try {
        triangleSet = new TriangleSet<T>(stlFileName, DBL);
    }
    catch (PlbIOException& exception) {
        pcout << "Error, could not read STL file " << stlFileName
              << ": " << exception.what() << std::endl;
        exit(-1);
    }

    if (readMaxIter) {
        triangleSet->refineByAreaRecursively(threshold, maxIter);
    } else {
        if (readThreshold) {
            triangleSet->refineByArea(threshold);
        } else {
            triangleSet->refine();
        }
    }

    try {
        triangleSet->writeAsciiSTL(outFileName);
    }
    catch (PlbIOException& exception) {
            pcout << "Error, could not write STL file " << outFileName
                  << ": " << exception.what() << std::endl;
            exit(-1);
    }

    delete triangleSet;

    return 0;
}

