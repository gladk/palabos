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

#ifndef TRIANGLE_SET_GENERATOR_HH
#define TRIANGLE_SET_GENERATOR_HH

#include "core/globalDefs.h"
#include "offLattice/triangleSetGenerator.h"

namespace plb {

template<typename T>
TriangleSet<T> constructSphere(Array<T,3> const& center, T radius, plint minNumOfTriangles)
{
    std::vector<typename TriangleSet<T>::Triangle> triangles;
#ifdef PLB_DEBUG
    static const T eps = std::numeric_limits<T>::epsilon();
#endif
    PLB_ASSERT(radius > (T) 0.0 && !util::fpequal(radius, (T) 0.0, eps) &&
               minNumOfTriangles >= 8);

    // Create a triangularized unit sphere

    // Initial 6 vertices

    Array<T,3> va;
    va[0] = (T) 1.0;
    va[1] = (T) 0.0;
    va[2] = (T) 0.0;

    Array<T,3> vb;
    vb[0] = (T) 0.0;
    vb[1] = (T) 1.0;
    vb[2] = (T) 0.0;

    Array<T,3> vc;
    vc[0] = (T) (-1.0);
    vc[1] = (T) 0.0;
    vc[2] = (T) 0.0;

    Array<T,3> vd;
    vd[0] = (T) 0.0;
    vd[1] = (T) (-1.0);
    vd[2] = (T) 0.0;

    Array<T,3> ve;
    ve[0] = (T) 0.0;
    ve[1] = (T) 0.0;
    ve[2] = (T) 1.0;

    Array<T,3> vf;
    vf[0] = (T) 0.0;
    vf[1] = (T) 0.0;
    vf[2] = (T) (-1.0);

    // Initial 8 triangles

    typename TriangleSet<T>::Triangle tmp;

    tmp[0] = ve;
    tmp[1] = va;
    tmp[2] = vb;
    triangles.push_back(tmp);

    tmp[0] = ve;
    tmp[1] = vb;
    tmp[2] = vc;
    triangles.push_back(tmp);

    tmp[0] = ve;
    tmp[1] = vc;
    tmp[2] = vd;
    triangles.push_back(tmp);

    tmp[0] = ve;
    tmp[1] = vd;
    tmp[2] = va;
    triangles.push_back(tmp);

    tmp[0] = vf;
    tmp[1] = vb;
    tmp[2] = va;
    triangles.push_back(tmp);

    tmp[0] = vf;
    tmp[1] = vc;
    tmp[2] = vb;
    triangles.push_back(tmp);

    tmp[0] = vf;
    tmp[1] = vd;
    tmp[2] = vc;
    triangles.push_back(tmp);

    tmp[0] = vf;
    tmp[1] = va;
    tmp[2] = vd;
    triangles.push_back(tmp);

    // Perform refinement iterations

    plint size;
    while ((size = triangles.size()) < minNumOfTriangles) {
        for (plint i = 0; i < size; i++) {
            va = triangles[i][0];
            vb = triangles[i][1];
            vc = triangles[i][2];

            vd = (T) 0.5 * (va + vb);
            ve = (T) 0.5 * (vb + vc);
            vf = (T) 0.5 * (vc + va);

            vd /= norm(vd);
            ve /= norm(ve);
            vf /= norm(vf);

            triangles[i][0] = vd;
            triangles[i][1] = ve;
            triangles[i][2] = vf;

            tmp[0] = va;
            tmp[1] = vd;
            tmp[2] = vf;
            triangles.push_back(tmp);

            tmp[0] = vd;
            tmp[1] = vb;
            tmp[2] = ve;
            triangles.push_back(tmp);

            tmp[0] = vf;
            tmp[1] = ve;
            tmp[2] = vc;
            triangles.push_back(tmp);
        }
    }

    // Scale and translate the mesh

    TriangleSet<T> triangleSet(triangles);

    triangleSet.scale(radius);
    triangleSet.translate(center);

    return triangleSet;
}

template<typename T>
TriangleSet<T> constructCylinder( Array<T,3> const& inletCenter, T inletRadius, T outletRadius,
                                  T length, plint nAxial, plint nCirc,
                                  std::vector<Array<T,3> >& inletPoints )
{
    static const T eps = std::numeric_limits<T>::epsilon();
    std::vector<typename TriangleSet<T>::Triangle> triangles;
    PLB_ASSERT( inletRadius > (T) 0.0 && !util::fpequal( inletRadius, (T) 0.0, eps) &&
               outletRadius > (T) 0.0 && !util::fpequal(outletRadius, (T) 0.0, eps) &&
                     length > (T) 0.0 && !util::fpequal(      length, (T) 0.0, eps) &&
               nAxial >= 2 && nCirc >= 2);

    inletPoints.resize(0);
    triangles.resize(0);

    // Construction of the cylindrical grid and triangulation

    T pi = std::acos((T) -1.0);

    T dtheta = (T) 2.0*pi / nCirc;
    T dx = length / (nAxial-(T)1.0);
    T dr = (outletRadius-inletRadius) / length * dx;

    T x0 = inletCenter[0];
    T y0 = inletCenter[1];
    T z0 = inletCenter[2];

    for (plint i = 0; i < nAxial-1; i++) {
        T x = x0 + i*dx;
        T r = (outletRadius-inletRadius) * (x-x0) / length + inletRadius;
        for (plint j = 0; j < nCirc; j++) {
            T theta = j * dtheta;

            Array<Array<T,3>, 4> v;

            v[0][0] = x;
            v[0][1] = y0 + r * std::cos(theta);
            v[0][2] = z0 + r * std::sin(theta);

            v[1][0] = x;
            v[1][1] = y0 + r * std::cos(theta + dtheta);
            v[1][2] = z0 + r * std::sin(theta + dtheta);

            v[2][0] = x + dx;
            v[2][1] = y0 + (r + dr) * std::cos(theta + dtheta);
            v[2][2] = z0 + (r + dr) * std::sin(theta + dtheta);

            v[3][0] = x + dx;
            v[3][1] = y0 + (r + dr) * std::cos(theta);
            v[3][2] = z0 + (r + dr) * std::sin(theta);

            Array<T,3> vc;

            vc[0] = x + (T) 0.5 * dx;
            vc[1] = y0 + (r + (T) 0.5 * dr) * std::cos(theta + (T) 0.5 * dtheta);
            vc[2] = z0 + (r + (T) 0.5 * dr) * std::sin(theta + (T) 0.5 * dtheta);

            Array<Array<T,3>, 4> vce;


            vce[0][0] = x;
            vce[0][1] = y0 + r * std::cos(theta + (T) 0.5 * dtheta);
            vce[0][2] = z0 + r * std::sin(theta + (T) 0.5 * dtheta);

            if (i==0) {
                inletPoints.push_back(v[0]);
                inletPoints.push_back(vce[0]);
            }

            vce[1][0] = x + (T) 0.5 * dx;
            vce[1][1] = y0 + (r + (T) 0.5 * dr) * std::cos(theta + dtheta);
            vce[1][2] = z0 + (r + (T) 0.5 * dr) * std::sin(theta + dtheta);

            vce[2][0] = x + dx;
            vce[2][1] = y0 + (r + dr) * std::cos(theta + (T) 0.5 * dtheta);
            vce[2][2] = z0 + (r + dr) * std::sin(theta + (T) 0.5 * dtheta);

            vce[3][0] = x + (T) 0.5 * dx;
            vce[3][1] = y0 + (r + (T) 0.5 * dr) * std::cos(theta);
            vce[3][2] = z0 + (r + (T) 0.5 * dr) * std::sin(theta);

            typename TriangleSet<T>::Triangle tmp;

            tmp[0] = vc;
            tmp[1] = v[0];
            tmp[2] = vce[0];
            triangles.push_back(tmp);

            tmp[0] = vc;
            tmp[1] = vce[0];
            tmp[2] = v[1];
            triangles.push_back(tmp);

            tmp[0] = vc;
            tmp[1] = v[1];
            tmp[2] = vce[1];
            triangles.push_back(tmp);

            tmp[0] = vc;
            tmp[1] = vce[1];
            tmp[2] = v[2];
            triangles.push_back(tmp);

            tmp[0] = vc;
            tmp[1] = v[2];
            tmp[2] = vce[2];
            triangles.push_back(tmp);

            tmp[0] = vc;
            tmp[1] = vce[2];
            tmp[2] = v[3];
            triangles.push_back(tmp);

            tmp[0] = vc;
            tmp[1] = v[3];
            tmp[2] = vce[3];
            triangles.push_back(tmp);

            tmp[0] = vc;
            tmp[1] = vce[3];
            tmp[2] = v[0];
            triangles.push_back(tmp);
        }
    }
    return TriangleSet<T>(triangles);
}
template<typename T>
TriangleSet<T> constructCylinder(Array<T,3> const& inletCenter, T inletRadius, T outletRadius,
                                 T length, plint nAxial, plint nCirc)
{
    static const T eps = std::numeric_limits<T>::epsilon();
    std::vector<typename TriangleSet<T>::Triangle> triangles;
    PLB_ASSERT( inletRadius > (T) 0.0 && !util::fpequal( inletRadius, (T) 0.0, eps) &&
               outletRadius > (T) 0.0 && !util::fpequal(outletRadius, (T) 0.0, eps) &&
                     length > (T) 0.0 && !util::fpequal(      length, (T) 0.0, eps) &&
               nAxial >= 2 && nCirc >= 2);

    triangles.resize(0);

    // Construction of the cylindrical grid and triangulation

    T pi = std::acos((T) -1.0);

    T dtheta = (T) 2.0*pi / nCirc;
    T dx = length / (nAxial-(T)1.0);
    T dr = (outletRadius-inletRadius) / length * dx;

    T x0 = inletCenter[0];
    T y0 = inletCenter[1];
    T z0 = inletCenter[2];

    for (plint i = 0; i < nAxial-1; i++) {
        T x = x0 + i*dx;
        T r = (outletRadius-inletRadius) * (x-x0) / length + inletRadius;
        for (plint j = 0; j < nCirc; j++) {
            T theta = j * dtheta;

            Array<Array<T,3>, 4> v;

            v[0][0] = x;
            v[0][1] = y0 + r * std::cos(theta);
            v[0][2] = z0 + r * std::sin(theta);

            v[1][0] = x;
            v[1][1] = y0 + r * std::cos(theta + dtheta);
            v[1][2] = z0 + r * std::sin(theta + dtheta);

            v[2][0] = x + dx;
            v[2][1] = y0 + (r + dr) * std::cos(theta + dtheta);
            v[2][2] = z0 + (r + dr) * std::sin(theta + dtheta);

            v[3][0] = x + dx;
            v[3][1] = y0 + (r + dr) * std::cos(theta);
            v[3][2] = z0 + (r + dr) * std::sin(theta);

            Array<T,3> vc;

            vc[0] = x + (T) 0.5 * dx;
            vc[1] = y0 + (r + (T) 0.5 * dr) * std::cos(theta + (T) 0.5 * dtheta);
            vc[2] = z0 + (r + (T) 0.5 * dr) * std::sin(theta + (T) 0.5 * dtheta);

            Array<Array<T,3>, 4> vce;

            vce[0][0] = x;
            vce[0][1] = y0 + r * std::cos(theta + (T) 0.5 * dtheta);
            vce[0][2] = z0 + r * std::sin(theta + (T) 0.5 * dtheta);

            vce[1][0] = x + (T) 0.5 * dx;
            vce[1][1] = y0 + (r + (T) 0.5 * dr) * std::cos(theta + dtheta);
            vce[1][2] = z0 + (r + (T) 0.5 * dr) * std::sin(theta + dtheta);

            vce[2][0] = x + dx;
            vce[2][1] = y0 + (r + dr) * std::cos(theta + (T) 0.5 * dtheta);
            vce[2][2] = z0 + (r + dr) * std::sin(theta + (T) 0.5 * dtheta);

            vce[3][0] = x + (T) 0.5 * dx;
            vce[3][1] = y0 + (r + (T) 0.5 * dr) * std::cos(theta);
            vce[3][2] = z0 + (r + (T) 0.5 * dr) * std::sin(theta);

            typename TriangleSet<T>::Triangle tmp;

            tmp[0] = vc;
            tmp[1] = v[0];
            tmp[2] = vce[0];
            triangles.push_back(tmp);

            tmp[0] = vc;
            tmp[1] = vce[0];
            tmp[2] = v[1];
            triangles.push_back(tmp);

            tmp[0] = vc;
            tmp[1] = v[1];
            tmp[2] = vce[1];
            triangles.push_back(tmp);

            tmp[0] = vc;
            tmp[1] = vce[1];
            tmp[2] = v[2];
            triangles.push_back(tmp);

            tmp[0] = vc;
            tmp[1] = v[2];
            tmp[2] = vce[2];
            triangles.push_back(tmp);

            tmp[0] = vc;
            tmp[1] = vce[2];
            tmp[2] = v[3];
            triangles.push_back(tmp);

            tmp[0] = vc;
            tmp[1] = v[3];
            tmp[2] = vce[3];
            triangles.push_back(tmp);

            tmp[0] = vc;
            tmp[1] = vce[3];
            tmp[2] = v[0];
            triangles.push_back(tmp);
        }
    }
    return TriangleSet<T>(triangles);
}

template<typename T>
void addSurface (
        Array<T,3> const& lowerCorner,
        Array<T,3> const& delta1, plint n1, Array<T,3> const& delta2, plint n2,
        std::vector<typename TriangleSet<T>::Triangle>& triangles )
{
    Array<T,3> pos1(lowerCorner);
    for(plint i1=0; i1<n1; ++i1, pos1+=delta1) {
        Array<T,3> pos2(pos1);
        for(plint i2=0; i2<n2; ++i2, pos2+=delta2) {
            typename TriangleSet<T>::Triangle triangle;
            triangle[0] = pos2;
            triangle[1] = pos2+delta1;
            triangle[2] = pos2+delta2;
            triangles.push_back(triangle);
            triangle[0] += delta1+delta2;
            std::swap(triangle[1], triangle[2]);
            triangles.push_back(triangle);
        }
    }
}

template<typename T>
TriangleSet<T> constructCuboid (
        Array<T,3> const& lowerCorner, Array<T,3> const& upperCorner,
        Array<plint,3> const& nSegments )
{
    std::vector<typename TriangleSet<T>::Triangle> triangles;
    T lx = upperCorner[0]-lowerCorner[0];
    T ly = upperCorner[1]-lowerCorner[1];
    T lz = upperCorner[2]-lowerCorner[2];
    Array<T,3> deltaX(lx/(T)nSegments[0], T(), T());
    Array<T,3> deltaY(T(), ly/(T)nSegments[1], T());
    Array<T,3> deltaZ(T(), T(), lz/(T)nSegments[2]);

    addSurface(lowerCorner,
               deltaZ, nSegments[2], deltaY, nSegments[1], triangles);
    addSurface(lowerCorner+Array<T,3>(lx,T(),T()),
               deltaY, nSegments[1], deltaZ, nSegments[2], triangles);

    addSurface(lowerCorner,
               deltaX, nSegments[0], deltaZ, nSegments[2], triangles);
    addSurface(lowerCorner+Array<T,3>(T(),ly,T()),
               deltaZ, nSegments[2], deltaX, nSegments[0], triangles);

    addSurface(lowerCorner,
               deltaY, nSegments[1], deltaX, nSegments[0], triangles);
    addSurface(lowerCorner+Array<T,3>(T(),T(),lz),
               deltaX, nSegments[0], deltaY, nSegments[1], triangles);

    return TriangleSet<T>(triangles);
}


template<typename T>
TriangleSet<T> constructCylinder( Array<T,3> const& inletCenter, Array<T,3> const& axis,
                                  T inletRadius, T outletRadius,
                                  T length, plint nAxial, plint nCirc,
                                  std::vector<Array<T,3> >& inletPoints )
{
    TriangleSet<T> cylinder (
            constructCylinder(Array<T,3>(T(),T(),T()), inletRadius, outletRadius,
            length, nAxial, nCirc, inletPoints) );
    Array<T,3> xAxis((T)1,T(),T());
    Array<T,3> rotAxis;
    T alignment = dot(xAxis, axis);
    T eps = (T)100.*std::numeric_limits<T>::epsilon();
    if (!util::fpequal(alignment, (T)1., eps)) {
        crossProduct<T>(xAxis, axis, rotAxis);
        rotAxis /= norm(rotAxis);
        T angle = angleBetweenVectors(xAxis, axis);
        cylinder.rotateAtOrigin(rotAxis, angle);
        for (pluint i=0; i<inletPoints.size(); ++i) {
            inletPoints[i] = rotateAtOrigin(inletPoints[i], rotAxis, angle);
        }
    }
    cylinder.translate(inletCenter);
    for (pluint i=0; i<inletPoints.size(); ++i) {
        inletPoints[i] += inletCenter;
    }
    return cylinder;
}


template<typename T>
TriangleSet<T> patchTubes(TriangleSet<T> const& geometryWithOpenings, plint sortDirection, std::vector<T> patchLengths)
{
    typedef typename TriangleSet<T>::Triangle Triangle;

    std::vector<Triangle> fullGeometry( geometryWithOpenings.getTriangles() );

    DEFscaledMesh<T>* defMesh = new DEFscaledMesh<T>(geometryWithOpenings);
    TriangularSurfaceMesh<T>& mesh = defMesh->getMesh();

    std::vector<Lid> holes = mesh.closeHoles();
    std::sort(holes.begin(), holes.end(), LidLessThan<T>(sortDirection, mesh));

    PLB_ASSERT( holes.size() == patchLengths.size() );
    
    for (pluint iHole=0; iHole<holes.size(); ++iHole) {
        Array<T,3> baryCenter = computeGeometricCenter(mesh,holes[iHole]);
        plint numHoleVertices = (plint) holes[iHole].boundaryVertices.size();

        Array<T,3> normal = computeNormal(mesh, holes[iHole]);
        T aveRadius = computeGeometricRadius(mesh,holes[iHole]);

        Array<T,3> nextCenter = baryCenter + normal*aveRadius*(T)10./(T)numHoleVertices;

        plint numInletPoints = numHoleVertices;
        bool oddNumber = numHoleVertices%2==1;
        if (oddNumber) numInletPoints--; // Must be even for cylinder construction algorithm.
        std::vector<Array<T,3> > inletPoints;
        plint numPointsOnLength = numInletPoints*patchLengths[iHole]/aveRadius/8;
        if (numPointsOnLength<3) numPointsOnLength = 3;
        TriangleSet<T> piece = constructCylinder(nextCenter, normal, aveRadius, aveRadius, patchLengths[iHole],
                                                 numPointsOnLength, numInletPoints/2, inletPoints);
        std::vector<Triangle> pieceTriangles = piece.getTriangles();

        plint newId = 0;
        T minDistance = std::numeric_limits<T>::max();
        plint minDistanceId = -1;
        for (plint i=0; i<numHoleVertices; ++i) {
            plint iVertex = holes[iHole].boundaryVertices[i];
            Array<T,3> p1 = mesh.getVertex(iVertex);
            T nextDistance = norm(inletPoints[newId]-p1);
            if (nextDistance<minDistance) {
                minDistance = nextDistance;
                minDistanceId = i;
            }
        }
        plint oldId = minDistanceId;
        plint newId_p1 = 0;
        for (plint i=0; i<numInletPoints; ++i) {
            plint newId_p1 = (newId+1) % numInletPoints;
            plint oldId_p1 = oldId-1;
            if (oldId_p1<0) oldId_p1 = numHoleVertices-1;

            plint oldVertex1 = holes[iHole].boundaryVertices[oldId];
            plint oldVertex2 = holes[iHole].boundaryVertices[oldId_p1];
            Array<T,3> p1 = mesh.getVertex(oldVertex1);
            Array<T,3> p2 = inletPoints[newId];
            Array<T,3> p3 = mesh.getVertex(oldVertex2);
            Array<T,3> p4 = inletPoints[newId_p1];

            pieceTriangles.push_back(Triangle(p1,p3,p2));
            pieceTriangles.push_back(Triangle(p2,p3,p4));

            std::swap(newId, newId_p1);
            std::swap(oldId, oldId_p1);
        }

        if (oddNumber) {
            plint id_a = holes[iHole].boundaryVertices[oldId];
            plint oldId_p2 = oldId-1;
            if (oldId_p2<0) oldId_p2 = numHoleVertices-1;
            plint id_b = holes[iHole].boundaryVertices[oldId_p2];
            plint id_c = newId_p1;
            Array<T,3> a = mesh.getVertex(id_a);
            Array<T,3> b = mesh.getVertex(id_b);
            Array<T,3> c = inletPoints[id_c];
            pieceTriangles.push_back(Triangle(a,b,c));
        }


        fullGeometry.insert(fullGeometry.end(), pieceTriangles.begin(), pieceTriangles.end());
    }

    return TriangleSet<T>(fullGeometry, geometryWithOpenings.getPrecision());
}

template<typename T>
TriangleSet<T> constructRectangle(T lx, T ly, plint nx, plint ny)
{
#ifdef PLB_DEBUG
    static const T eps = std::numeric_limits<T>::epsilon();
#endif
    PLB_ASSERT(lx > (T) 0.0 && !util::fpequal(lx, (T) 0.0, eps) &&
               ly > (T) 0.0 && !util::fpequal(ly, (T) 0.0, eps) &&
               nx >= 2 && ny >= 2);

    T dx = lx / (T) (nx - 1);
    T dy = ly / (T) (ny - 1);

    std::vector<typename TriangleSet<T>::Triangle> triangles;
    T z = (T) 0;
    for (plint i = 0; i < nx-1; i++) {
        T x0 = i * dx;
        T x1 = x0 + dx;
        for (plint j = 0; j < ny-1; j++) {
            T y0 = j * dy;
            T y1 = y0 + dy;

            Array<T,3> v0, v1, v2;
            v0 = Array<T,3>(x0, y0, z);
            v1 = Array<T,3>(x1, y0, z);
            v2 = Array<T,3>(x1, y1, z);

            typename TriangleSet<T>::Triangle tmp;
            tmp[0] = v0;
            tmp[1] = v1;
            tmp[2] = v2;

            triangles.push_back(tmp);

            v1 = Array<T,3>(x0, y1, z);

            tmp[1] = v2;
            tmp[2] = v1;

            triangles.push_back(tmp);
        }
    }

    return TriangleSet<T>(triangles);
}

template<typename T>
TriangleSet<T> constructGenericRectangle(T lx, T ly, plint nx, plint ny, Array<T,3> const& center, Array<T,3> const& normal)
{
    static const T eps = std::numeric_limits<T>::epsilon();
    PLB_ASSERT(lx > (T) 0.0 && !util::fpequal(lx, (T) 0.0, eps) &&
               ly > (T) 0.0 && !util::fpequal(ly, (T) 0.0, eps) &&
               nx >= 2 && ny >= 2);

    TriangleSet<T> rectangle = constructRectangle<T>(lx, ly, nx, ny);
    Array<T,3> oldCenter((T)0.5 * lx, (T)0.5 * ly, (T) 0.0);
    rectangle.translate(-oldCenter);

    T r = norm(normal);
    PLB_ASSERT(!util::fpequal(r, (T) 0.0, eps));
    Array<T,3> unitNormal = normal / r;
    T theta = std::acos(unitNormal[2]);
    Array<T,3> zAxis((T) 0.0, (T) 0.0, (T) 1.0);
    Array<T,3> normedAxis = crossProduct<T>(zAxis, unitNormal);
    T l = norm(normedAxis);
    if (util::fpequal_abs(l, (T) 0.0, eps)) {
        normedAxis = Array<T,3>((T) 1.0, (T) 0.0, (T) 0.0);
        if (dot<T>(zAxis, unitNormal) > 0.0) {
            theta = 0.0;
        } else {
            theta = std::acos((T) -1.0);
        }
    } else {
        normedAxis /= l;
    }
    rectangle.rotateAtOrigin(normedAxis, theta);
    rectangle.translate(center);

    return rectangle;
}

template<typename T>
TriangleSet<T> constructStrip(std::vector<Array<T,3> > from, std::vector<Array<T,3> > to)
{
    PLB_ASSERT(from.size() >= 2);
    PLB_ASSERT(from.size() == to.size());
    plint size = from.size();

    std::vector<typename TriangleSet<T>::Triangle> triangles;

    plint count = 0;
    for (plint i = 0; i < size-1; i++) {
        Array<T,3> from0 = from[i];
        Array<T,3> from1 = from[i+1];

        Array<T,3> to0 = to[i];
        Array<T,3> to1 = to[i+1];

        typename TriangleSet<T>::Triangle triangle;
        if (count % 2 == 0) {
            triangle[0] = from0;
            triangle[1] = from1;
            triangle[2] = to0;
            triangles.push_back(triangle);

            triangle[0] = from1;
            triangle[1] = to1;
            triangle[2] = to0;
            triangles.push_back(triangle);
        } else {
            triangle[0] = from0;
            triangle[1] = from1;
            triangle[2] = to1;
            triangles.push_back(triangle);

            triangle[0] = from0;
            triangle[1] = to1;
            triangle[2] = to0;
            triangles.push_back(triangle);
        }
        count++;
    }

    return TriangleSet<T>(triangles);
}

} // namespace plb

#endif  // TRIANGLE_SET_GENERATOR_HH
