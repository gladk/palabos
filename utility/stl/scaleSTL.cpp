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
    plint fitDirection;
    T fitLength;
    try {
        global::argv(1).read(stlFileName);
        global::argv(2).read(fitDirection);
        global::argv(3).read(fitLength);
        global::argv(4).read(outFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " inputSTL.stl fitDirection fitLength outputSTL.stl" << std::endl;
        pcout << "If for example, you want your STL to have a length of 8m in x-direction (at the largest diameter), you should write" << std::endl;
        pcout << (std::string)global::argv(0) << " inputSTL.stl 0 8.0 outputSTL.stl" << std::endl;
        exit(-1);
    }

    if(fitDirection<0 || fitDirection>2) {
        pcout << "Error, fitDirection must be 0 (for x), 1 (for y), or 2 (for z)" << std::endl;
        exit(-1);
    }

    TriangleSet<T>* triangleSet = 0;
    try {
        triangleSet = new TriangleSet<T>(stlFileName, DBL);
    }
    catch (PlbIOException& exception) {
        pcout << "Error, could not read STL file " << stlFileName
              << ": " << exception.what() << std::endl;
        return -1;
    }
    try {
        DEFscaledMesh<T>* defMesh = new DEFscaledMesh<T>(*triangleSet);
        Array<T,2> xRange, yRange, zRange;
        defMesh->getMesh().computeBoundingBox(xRange, yRange, zRange);
        T scale = 0.;
        switch(fitDirection) {
            case 0:
                scale = fitLength / (xRange[1]-xRange[0]);
                break;
            case 1:
                scale = fitLength / (yRange[1]-yRange[0]);
                break;
            case 2:
                scale = fitLength / (zRange[1]-zRange[0]);
                break;
            default:
                PLB_ASSERT( false );
        }
        pcout << "The scale factor is " << scale << std::endl;
        defMesh->getMesh().scale(scale);
        defMesh->getMesh().writeAsciiSTL(outFileName);
    }
    catch (PlbIOException& exception) {
            pcout << "Error in STL file " << stlFileName
                  << ": " << exception.what() << std::endl;
            exit(-1);
    }

    return 0;
}

