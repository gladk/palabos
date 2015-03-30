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

typedef TriangleSet<T>::Triangle Triangle;

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./");
    global::IOpolicy().activateParallelIO(false);

    std::string precisionStr;
    try {
        global::argv(1).read(precisionStr);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " [FLT | DBL | LDBL | INF] inputSTL.stl" << std::endl;
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
    pcout << "STL file to TriangleSet SUCCESSFUL!" << std::endl;

    if (set->hasDegenerateTriangles()) {
        pcout << "WARNING: the TriangleSet has degenerate triangles!" << std::endl;
    } else {
        pcout << "The TriangleSet does not have degenerate triangles." << std::endl;
    }

    ConnectedTriangleSet<T>* connectedSet = 0;
    try {
        connectedSet = new ConnectedTriangleSet<T>(*set);
    }
    catch (PlbIOException& exception) {
        pcout << "ERROR, could not create a ConnectedTriangleSet object: "
              << exception.what() << std::endl;
        exit(-1);
    }
    pcout << "TriangleSet to ConnectedTriangleSet SUCCESSFUL!" << std::endl;
    pcout << "The connected mesh has: " << connectedSet->getNumVertices() << " vertices, and "
          << connectedSet->getNumTriangles() << " triangles." << std::endl;
    std::string offFileName(stlFileName);
    std::string extension = ".stl";
    size_t found = offFileName.rfind(extension);
    if (found != std::string::npos) {
        offFileName.replace(found, extension.length(), ".off");
        connectedSet->writeOFF(offFileName);
    }
    for (plint iTriangle = 0; iTriangle < connectedSet->getNumTriangles(); iTriangle++) {
        T area;
        Array<T,3> unitNormal;
        connectedSet->computeTriangleAreaAndUnitNormal(iTriangle, area, unitNormal);
        if (norm<T,3>(unitNormal) < (T) 0.5) {
            Array<plint,3> t = connectedSet->getTriangle(iTriangle);
            Array<T,3> v0 = connectedSet->getVertex(t[0]);
            Array<T,3> v1 = connectedSet->getVertex(t[1]);
            Array<T,3> v2 = connectedSet->getVertex(t[2]);
            pcout << "Computation of the unit triangle normal failed for the triangle with vertices: " << std::endl;
            pcout << "(" << v0[0] << ", " << v0[1] << ", " << v0[2] << ")" << std::endl;
            pcout << "(" << v1[0] << ", " << v1[1] << ", " << v1[2] << ")" << std::endl;
            pcout << "(" << v2[0] << ", " << v2[1] << ", " << v2[2] << ")" << std::endl;

            std::vector<Triangle> badTriangles;
            Triangle triangle;
            triangle[0] = v0;
            triangle[1] = v1;
            triangle[2] = v2;
            badTriangles.push_back(triangle);
            TriangleSet<T> badTriangleSet(badTriangles, precision);
            badTriangleSet.writeAsciiSTL("badTriangles.stl");
            exit(1);
        }
    }
    pcout << "Computation of triangle areas and unit normals SUCCESSFUL!" << std::endl;
    for (plint iVertex = 0; iVertex < connectedSet->getNumVertices(); iVertex++) {
        T area;
        Array<T,3> unitNormal;
        connectedSet->computeVertexAreaAndUnitNormal(iVertex, area, unitNormal);
        if (norm<T,3>(unitNormal) < (T) 0.5) {
            Array<T,3> v = connectedSet->getVertex(iVertex);
            pcout << "Computation of the unit vertex normal failed for the vertex: " << std::endl;
            pcout << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")" << std::endl;

            std::vector<plint> trianglesOnVertex = connectedSet->getTrianglesOnVertex(iVertex);
            std::vector<Triangle> badTriangles;
            for (plint i = 0; i < (plint) trianglesOnVertex.size(); i++) {
                plint iTriangle = trianglesOnVertex[i];
                Array<plint,3> t = connectedSet->getTriangle(iTriangle);
                Triangle triangle;
                triangle[0] = connectedSet->getVertex(t[0]);
                triangle[1] = connectedSet->getVertex(t[1]);
                triangle[2] = connectedSet->getVertex(t[2]);
                badTriangles.push_back(triangle);
            }
            TriangleSet<T> badTriangleSet(badTriangles, precision);
            badTriangleSet.writeAsciiSTL("badTriangles.stl");
            exit(1);
        }
    }
    pcout << "Computation of vertex areas and unit normals SUCCESSFUL!" << std::endl;

    DEFscaledMesh<T>* mesh = 0;
    try {
        mesh = new DEFscaledMesh<T>(*set);
    }
    catch (PlbIOException& exception) {
        pcout << "ERROR, could not create a DEFscaledMesh object: "
              << exception.what() << std::endl;
        exit(-1);
    }
    pcout << "TriangleSet to DEFscaledMesh SUCCESSFUL!" << std::endl;
    pcout << "The DEF mesh has: " << mesh->getMesh().getNumVertices() << " vertices, "
          << mesh->getMesh().getNumTriangles() << " triangles, and "
          << mesh->getMesh().detectHoles().size() << " holes." << std::endl;

    TriangleBoundary3D<T>* boundary = 0;
    try {
        boundary = new TriangleBoundary3D<T>(*mesh);
    }
    catch (PlbIOException& exception) {
        pcout << "ERROR, could not create a TriangleBoundary3D object: "
              << exception.what() << std::endl;
        exit(-1);
    }
    pcout << "DEFscaledMesh to TriangleBoundary3D SUCCESSFUL!" << std::endl;
    pcout << "The TriangleBoundary3D mesh has: " << boundary->getMesh().getNumVertices() << " vertices, and "
          << boundary->getMesh().getNumTriangles() << " triangles." << std::endl;

    pcout << std::endl << "All checks SUCCESSFUL!" << std::endl << std::endl;

    delete boundary;
    delete mesh;
    delete connectedSet;
    delete set;

    return 0;
}

