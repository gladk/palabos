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
#include <vector>

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
    std::string baseFileName(FileName(stlFileName).getName());

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

    int numWarnings = 0;

    plint numZeroAreaTriangles = set->numZeroAreaTriangles();
    if (numZeroAreaTriangles) {
        pcout << std::endl;
        pcout << "WARNING: Number of zero-area triangles in the TriangleSet: " << numZeroAreaTriangles << std::endl;
        pcout << std::endl;
        numWarnings++;
    } else {
        pcout << "The TriangleSet does not have any zero-area triangles." << std::endl;
    }

    if (set->hasFloatingPointPrecisionDependence()) {
        pcout << std::endl;
        pcout << "WARNING: The TriangleSet has a dependence on the floating point precision used." << std::endl;
        pcout << std::endl;
        numWarnings++;
    } else {
        pcout << "The TriangleSet does not have a dependence on the floating point precision used." << std::endl;
    }
    plint numTriangleSetTriangles = set->getTriangles().size();
    pcout << "TriangleSet:" << std::endl;
    pcout << "    Number of triangles: " << numTriangleSetTriangles << std::endl;

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
    plint numConnectedSetVertices = connectedSet->getNumVertices();
    plint numConnectedSetTriangles = connectedSet->getNumTriangles();
    pcout << "ConnectedTriangleSet:" << std::endl;
    pcout << "  Number of vertices : " << numConnectedSetVertices << std::endl;
    pcout << "  Number of triangles: " << numConnectedSetTriangles << std::endl;
    if (numConnectedSetTriangles != numTriangleSetTriangles) {
        pcout << std::endl;
        pcout << "WARNING: The TriangleSet and the ConnectedTriangleSet do not have the same number of triangles." << std::endl;
        pcout << std::endl;
        numWarnings++;
    }

    std::vector<Triangle> badTriangles;
    bool succeeded = true;
    for (plint iTriangle = 0; iTriangle < connectedSet->getNumTriangles(); iTriangle++) {
        T area;
        Array<T,3> unitNormal;
        connectedSet->computeTriangleAreaAndUnitNormal(iTriangle, area, unitNormal);
        if (norm<T,3>(unitNormal) < (T) 0.5) {
            succeeded = false;

            Array<plint,3> t = connectedSet->getTriangle(iTriangle);
            Array<T,3> v0 = connectedSet->getVertex(t[0]);
            Array<T,3> v1 = connectedSet->getVertex(t[1]);
            Array<T,3> v2 = connectedSet->getVertex(t[2]);
            pcout << "Computation of the unit triangle normal failed for the triangle with vertices: " << std::endl;
            pcout << "(" << v0[0] << ", " << v0[1] << ", " << v0[2] << ")" << std::endl;
            pcout << "(" << v1[0] << ", " << v1[1] << ", " << v1[2] << ")" << std::endl;
            pcout << "(" << v2[0] << ", " << v2[1] << ", " << v2[2] << ")" << std::endl;

            Triangle triangle;
            triangle[0] = v0;
            triangle[1] = v1;
            triangle[2] = v2;
            badTriangles.push_back(triangle);
        }
    }

    if (!succeeded) {
        pcout << std::endl;
        pcout << "WARNING: Number of triangles for which the computation of unit normals failed: " <<
           badTriangles.size() << std::endl;
        pcout << "         The problematic triangles are saved in the file badTriangles.stl" << std::endl;
        pcout << std::endl;
        numWarnings++;

        TriangleSet<T> badTriangleSet(badTriangles, precision);
        badTriangleSet.writeAsciiSTL("badTriangles.stl");
    } else {
        pcout << "Computation of triangle areas and unit normals SUCCESSFUL!" << std::endl;
    }

    badTriangles.clear();
    succeeded = true;
    plint numBadVertices = 0;
    for (plint iVertex = 0; iVertex < connectedSet->getNumVertices(); iVertex++) {
        T area;
        Array<T,3> unitNormal;
        connectedSet->computeVertexAreaAndUnitNormal(iVertex, area, unitNormal);
        if (norm<T,3>(unitNormal) < (T) 0.5) {
            succeeded = false;
            numBadVertices++;

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
        }
    }

    if (!succeeded) {
        pcout << std::endl;
        pcout << "WARNING: Number of veritces for which the computation of unit normals failed: " <<
           numBadVertices << std::endl;
        pcout << "         The associated triangles are saved in the file badVertices.stl" << std::endl;
        pcout << std::endl;
        numWarnings++;

        TriangleSet<T> badTriangleSet(badTriangles, precision);
        badTriangleSet.writeAsciiSTL("badVertices.stl");
    } else {
        pcout << "Computation of vertex areas and unit normals SUCCESSFUL!" << std::endl;
    }

    connectedSet->writeOFF(baseFileName + ".off");

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
    plint numDefVertices = mesh->getMesh().getNumVertices();
    plint numDefTriangles = mesh->getMesh().getNumTriangles();
    plint numHoles = mesh->getMesh().detectHoles().size();
    pcout << "DEFscaledMesh:" << std::endl;
    pcout << "  Number of vertices : " << numDefVertices << std::endl;
    pcout << "  Number of triangles: " << numDefTriangles << std::endl;
    pcout << "  Number of holes    : " << numHoles << std::endl;
    if (numDefVertices != numConnectedSetVertices) {
        pcout << std::endl;
        pcout << "WARNING: The ConnectedTriangleSet and the DEFscaledMesh do not have the same number of vertices." << std::endl;
        pcout << std::endl;
        numWarnings++;
    }
    if (numDefTriangles != numConnectedSetTriangles) {
        pcout << std::endl;
        pcout << "WARNING: The ConnectedTriangleSet and the DEFscaledMesh do not have the same number of triangles." << std::endl;
        pcout << std::endl;
        numWarnings++;
    }

    TriangleBoundary3D<T>* boundary = 0;
    try {
        bool automaticCloseHoles = false;
        boundary = new TriangleBoundary3D<T>(*mesh, automaticCloseHoles);
    }
    catch (PlbIOException& exception) {
        pcout << "ERROR, could not create a TriangleBoundary3D object: "
              << exception.what() << std::endl;
        exit(-1);
    }
    pcout << "DEFscaledMesh to TriangleBoundary3D SUCCESSFUL!" << std::endl;
    plint numTriangleBoundaryVertices = boundary->getMesh().getNumVertices();
    plint numTriangleBoundaryTriangles = boundary->getMesh().getNumTriangles();
    pcout << "TriangleBoundary3D:" << std::endl;
    pcout << "  Number of vertices : " << numTriangleBoundaryVertices << std::endl;
    pcout << "  Number of triangles: " << numTriangleBoundaryTriangles << std::endl;
    if (numTriangleBoundaryVertices != numConnectedSetVertices) {
        pcout << std::endl;
        pcout << "WARNING: The ConnectedTriangleSet and the TriangleBoundary3D do not have the same number of vertices." << std::endl;
        pcout << std::endl;
        numWarnings++;
    }
    if (numTriangleBoundaryTriangles != numConnectedSetTriangles) {
        pcout << std::endl;
        pcout << "WARNING: The ConnectedTriangleSet and the TriangleBoundary3D do not have the same number of triangles." << std::endl;
        pcout << std::endl;
        numWarnings++;
    }

    boundary->getMesh().writeAsciiSTL(baseFileName + "_TriangleBoundary3D.stl");

    TriangleBoundary3D<T>* closedBoundary = 0;
    if (numHoles) {
        try {
            bool automaticCloseHoles = true;
            closedBoundary = new TriangleBoundary3D<T>(*mesh, automaticCloseHoles);
        }
        catch (PlbIOException& exception) {
            pcout << "ERROR, could not create a closed TriangleBoundary3D object: "
                  << exception.what() << std::endl;
            exit(-1);
        }
        pcout << "DEFscaledMesh to a closed TriangleBoundary3D SUCCESSFUL!" << std::endl;
        plint numClosedTriangleBoundaryVertices = closedBoundary->getMesh().getNumVertices();
        plint numClosedTriangleBoundaryTriangles = closedBoundary->getMesh().getNumTriangles();
        pcout << "Closed TriangleBoundary3D:" << std::endl;
        pcout << "  Number of vertices : " << numClosedTriangleBoundaryVertices << std::endl;
        pcout << "  Number of triangles: " << numClosedTriangleBoundaryTriangles << std::endl;

        closedBoundary->getMesh().writeAsciiSTL(baseFileName + "_closed_TriangleBoundary3D.stl");
    }

    pcout << std::endl << "All checks SUCCESSFUL! ";
    if (numWarnings) {
        pcout << "Number of WARNINGS issued: " << numWarnings;
    }
    pcout << std::endl << std::endl;

    delete closedBoundary;
    delete boundary;
    delete mesh;
    delete connectedSet;
    delete set;

    return 0;
}
