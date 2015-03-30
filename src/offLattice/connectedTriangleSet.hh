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

/* Main author: Dimitrios Kontaxakis */

#ifndef CONNECTED_TRIANGLE_SET_HH
#define CONNECTED_TRIANGLE_SET_HH

#include "core/array.h"
#include "core/globalDefs.h"
#include "offLattice/triangleSet.h"
#include "offLattice/connectedTriangleSet.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

#include <cmath>
#include <cstdio>

namespace plb {

/* ******** Class ConnectedTriangleSet::VertexLessThan *************************** */

template<typename T>
inline bool ConnectedTriangleSet<T>::VertexSetLessThan::vertexComponentLessThan(T x, T y)
{
    T tmp = (T) 0.5 * (std::fabs(x - y) + std::fabs(y - x));
    return ((x < y) && (tmp > epsilon));
}

template<typename T>
inline bool ConnectedTriangleSet<T>::VertexSetLessThan::vertexComponentEqual(T x, T y)
{
    return ((!vertexComponentLessThan(x, y)) && (!vertexComponentLessThan(y, x)));
}

template<typename T>
inline bool ConnectedTriangleSet<T>::VertexSetLessThan::vertexLessThan(
        Array<T,3> const& v1, Array<T,3> const& v2)
{
    return ((vertexComponentLessThan(v1[0], v2[0]) ||
            (vertexComponentEqual(v1[0], v2[0]) && vertexComponentLessThan(v1[1], v2[1])) ||
            (vertexComponentEqual(v1[0], v2[0]) && vertexComponentEqual(v1[1], v2[1])
                                                && vertexComponentLessThan(v1[2], v2[2]))));
}

template<typename T>
inline bool ConnectedTriangleSet<T>::VertexSetLessThan::operator()(
        VertexSetNode const& node1, VertexSetNode const& node2)
{
    Array<T,3> const& v1 = *(node1.vertex);
    Array<T,3> const& v2 = *(node2.vertex);

    return vertexLessThan(v1, v2);
}

/* ******** Class ConnectedTriangleSet ******************************************* */

template<typename T>
ConnectedTriangleSet<T>::ConnectedTriangleSet(TriangleSet<T> const& triangleSet)
{
    std::vector<Array<Array<T,3>,3> > oldTriangles = triangleSet.getTriangles();

    numTriangles = (plint) oldTriangles.size();
    if (numTriangles == 0) {
        return;
    }

    triangles.resize(numTriangles);

    T epsilon = getEpsilon<T>(triangleSet.getPrecision());
    VertexSetLessThan lessThan(epsilon);
    VertexSet vertexSet(lessThan);

    // Unify duplicated vertices.
    numVertices = 0;
    plint globalVertex = 0;
    for (plint iTriangle = 0; iTriangle < numTriangles; iTriangle++) {
        for (plint localVertex = 0; localVertex < 3; localVertex++) {
            Array<T,3> *vertex = &oldTriangles[iTriangle][localVertex];
            VertexSetIterator it = vertexSet.find(VertexSetNode(-1, vertex));
            if (it == vertexSet.end()) {
                globalVertex = numVertices;
                numVertices++;
                VertexSetNode newNode(globalVertex, vertex);
                vertexSet.insert(newNode);
            } else {
                globalVertex = it->i;
            }
            triangles[iTriangle][localVertex] = globalVertex;
        }
    }
    PLB_ASSERT(numVertices >= 3);

    // Order and copy vertex coordinates.
    vertices.resize(numVertices);
    VertexSetConstIterator it = vertexSet.begin();
    for (; it != vertexSet.end(); ++it) {
        vertices[it->i] = *(it->vertex);
    }

    // Create the connectivity information on the triangle global indices
    // that meet on a unified vertex.
    trianglesOnVertex.resize(numVertices);
    for (plint iTriangle = 0; iTriangle < numTriangles; iTriangle++) {
        for (plint localVertex = 0; localVertex < 3; localVertex++) {
            plint globalVertex = triangles[iTriangle][localVertex];
            trianglesOnVertex[globalVertex].push_back(iTriangle);
        }
    }
}

template<typename T>
void ConnectedTriangleSet<T>::swapGeometry(std::vector<Array<T,3> >& newVertices)
{
    PLB_ASSERT(newVertices.size() == vertices.size());
    vertices.swap(newVertices);
}

template<typename T>
void ConnectedTriangleSet<T>::computeVertexAreaAndUnitNormal(plint iVertex, T& area,
        Array<T,3>& unitNormal, std::vector<Array<T,3> > *newVertices, plint indexOffset) const
{
    PLB_ASSERT(iVertex >= 0 && iVertex < numVertices);
    area = (T) 0;
    unitNormal.resetToZero();
    for (plint iTriangle = 0; iTriangle < (plint) trianglesOnVertex[iVertex].size(); iTriangle++) {
        plint globalTriangle = trianglesOnVertex[iVertex][iTriangle];
        T triangleArea = (T) 0;
        Array<T,3> triangleNormal((T) 0, (T) 0, (T) 0);
        computeTriangleAreaAndUnitNormal(globalTriangle, triangleArea, triangleNormal, newVertices, indexOffset);
        area += triangleArea;
        unitNormal += triangleNormal;
    }
    area /= (T) 3;

    T normNormal = norm<T,3>(unitNormal);
    if (normNormal <= getEpsilon<T>(floatingPointPrecision<T>())) {
        unitNormal.resetToZero();
    } else {
        unitNormal /= normNormal;
    }
}

template<typename T>
void ConnectedTriangleSet<T>::computeTriangleAreaAndUnitNormal(plint iTriangle, T& area,
        Array<T,3>& unitNormal, std::vector<Array<T,3> > *newVertices, plint indexOffset) const
{
    PLB_ASSERT(iTriangle >= 0 && iTriangle < numTriangles);

    std::vector<Array<T,3> > const* workingVertices = newVertices;
    if (workingVertices == 0) {
        workingVertices = &vertices;
        indexOffset = 0;
    }

    plint v0 = triangles[iTriangle][0] + indexOffset;
    plint v1 = triangles[iTriangle][1] + indexOffset;
    plint v2 = triangles[iTriangle][2] + indexOffset;
    Array<T,3> edge01 = (*workingVertices)[v1] - (*workingVertices)[v0];
    Array<T,3> edge02 = (*workingVertices)[v2] - (*workingVertices)[v0];
    crossProduct<T>(edge01, edge02, unitNormal);
    T normNormal = norm<T,3>(unitNormal);
    if (normNormal <= getEpsilon<T>(floatingPointPrecision<T>())) {
        unitNormal.resetToZero();
        area = (T) 0;
        return;
    }
    unitNormal /= normNormal;
    area = 0.5 * normNormal;
}

template<typename T>
void ConnectedTriangleSet<T>::writeOFF(std::string fname, std::vector<Array<T,3> > *newVertices,
        plint indexOffset) const
{
    if (global::mpi().isMainProcessor()) {
        if (numTriangles == 0) {
            return;
        }
        FILE *fp = fopen(fname.c_str(), "w");
        PLB_ASSERT(fp != NULL);

        std::vector<Array<T,3> > const* workingVertices = newVertices;
        if (workingVertices == 0) {
            workingVertices = &vertices;
            indexOffset = 0;
        }

        fprintf(fp, "OFF\n");
        fprintf(fp, "%ld %ld %ld\n", (long) numVertices, (long) numTriangles, (long) 0);
        for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
            fprintf(fp, "% .10e % .10e % .10e\n", (double) (*workingVertices)[iVertex+indexOffset][0],
                    (double) (*workingVertices)[iVertex+indexOffset][1],
                    (double) (*workingVertices)[iVertex+indexOffset][2]);
        }
        for (plint iTriangle = 0; iTriangle < numTriangles; iTriangle++) {
            fprintf(fp, "%ld %ld %ld %ld\n", (long) 3, (long) triangles[iTriangle][0],
                    (long) triangles[iTriangle][1],  (long) triangles[iTriangle][2]);
        }
    }
}

template<typename T>
TriangleSet<T>* ConnectedTriangleSet<T>::toTriangleSet(Precision precision, std::vector<Array<T,3> > *newVertices,
        plint indexOffset) const
{
    if (numTriangles == 0) {
        return 0;
    }

    std::vector<Array<T,3> > const* workingVertices = newVertices;
    if (workingVertices == 0) {
        workingVertices = &vertices;
        indexOffset = 0;
    }

    typedef typename TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> vectorOfTriangles;

    for (plint iTriangle = 0; iTriangle < numTriangles; iTriangle++) {
        plint iVertex0 = triangles[iTriangle][0] + indexOffset;
        plint iVertex1 = triangles[iTriangle][1] + indexOffset;
        plint iVertex2 = triangles[iTriangle][2] + indexOffset;
        Array<T,3> vertex0 = (*workingVertices)[iVertex0];
        Array<T,3> vertex1 = (*workingVertices)[iVertex1];
        Array<T,3> vertex2 = (*workingVertices)[iVertex2];
        Triangle triangle(vertex0, vertex1, vertex2);
        vectorOfTriangles.push_back(triangle);
    }

    return new TriangleSet<T>(vectorOfTriangles, precision);
}

} // namespace plb

#endif  // CONNECTED_TRIANGLE_SET_HH
