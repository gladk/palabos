/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2012 FlowKit Sarl
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
#include "palabos2D.h"
#include "palabos2D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor


template<typename T>
class ShapeFunction {
public:
    ShapeFunction(Array<plint,3> const& c1_, Array<plint,3> const& c2_)
        : c1(c1_), c2(c2_)
    { }
    T operator()(plint iX, plint iY, plint iZ) {
        return std::min( sqrt(util::sqr(iX-c1[0])+util::sqr(iY-c1[1])+util::sqr(iZ-c1[2])),
                         sqrt(util::sqr(iX-c2[0])+util::sqr(iY-c2[1])+util::sqr(iZ-c2[2])) );
    }
private:
    Array<plint,3> c1, c2;
};


int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);

    // First application: Calculate two iso-surfaces from a scalar-field.
    plint nx=200, ny=200, nz=200;
    std::vector<T> isoLevels;
    isoLevels.push_back(50);
    isoLevels.push_back(100);
    MultiScalarField3D<T> scalarField(nx,ny,nz);
    setToFunction(scalarField, scalarField.getBoundingBox(), ShapeFunction<T>(Array<plint,3>(40,40,40),Array<plint,3>(140,140,140)));

    typedef TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> triangles;
    isoSurfaceMarchingCube(triangles, scalarField, isoLevels, scalarField.getBoundingBox().enlarge(-1));
    TriangleSet<T>(triangles).writeAsciiSTL("iso.stl");

    // Second application: Revoxelize an STL file.
    TriangleSet<T> artery("aneurysm.stl", DBL);
    plint blockSize = 20;
    plint resolution = 80;
    plint referenceDirection = 0;

    plint margin = 1;
    plint extraLayer = 0;
    plint borderWidth = 1;
    plint envelopeWidth = 1;
    DEFscaledMesh<T>* defMesh =
                new DEFscaledMesh<T>(artery, resolution, referenceDirection, margin, extraLayer);
    TriangleBoundary3D<T> boundary(*defMesh);
    delete defMesh;
    boundary.getMesh().inflate();
    const int flowType = voxelFlag::inside;
    VoxelizedDomain3D<T> voxelizedDomain (
                boundary, flowType, extraLayer, borderWidth, envelopeWidth, blockSize );

    isoSurfaceMarchingCube(triangles, voxelizedDomain, voxelizedDomain.getVoxelMatrix().getBoundingBox());

    // Transform the created surface back from grid coordinates to original coordinates.
    TriangleSet<T> triangleSet(triangles);
    triangleSet.scale(voxelizedDomain.getBoundary().getDx());
    triangleSet.translate(voxelizedDomain.getBoundary().getPhysicalLocation());
    triangleSet.writeAsciiSTL("newartery.stl");
}

