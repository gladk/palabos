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
    Array<T,3> normedAxis;
    T theta;
    try {
        global::argv(1).read(stlFileName);
        global::argv(2).read(normedAxis[0]);
        global::argv(3).read(normedAxis[1]);
        global::argv(4).read(normedAxis[2]);
        global::argv(5).read(theta);
        global::argv(6).read(outFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " inputSTL.stl normedAxis_x normedAxis_y normedAxis_z angle outputSTL.stl"
              << std::endl;
        pcout << "If for example, you want your STL to be rotated around the y-axis for an angle of 12 degrees, you should write"
              << std::endl;
        pcout << (std::string)global::argv(0) << " inputSTL.stl 0.0 1.0 0.0 12.0 outputSTL.stl" << std::endl;
        exit(-1);
    }

    T normNormedAxis = norm(normedAxis);
    if (util::isZero(normNormedAxis)) {
        pcout << "The normedAxis should be a vector of magnitude 1." << std::endl;
        return -1;
    }
    normedAxis /= norm(normedAxis);

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
        T pi = std::acos((T) -1);
        theta *= (pi/(T)180.0);
        triangleSet->rotateAtOrigin(normedAxis, theta);
        triangleSet->writeAsciiSTL(outFileName);
    }
    catch (PlbIOException& exception) {
            pcout << "Error in STL file " << stlFileName
                  << ": " << exception.what() << std::endl;
            exit(-1);
    }

    return 0;
}

