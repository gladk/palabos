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

#include "parallelism/mpiManager.h"
#include "core/globalDefs.h"
#include "core/serializer.h"
#include "core/serializer.hh"
#include "io/vtkDataOutput.h"
#include "io/vtkDataOutput.hh"
#include "io/serializerIO.h"
#include "io/multiBlockWriter3D.h"
#include "io/base64.h"
#include "io/base64.hh"
#include <fstream>

namespace plb {
    
////////// class VtkDataWriter3D ////////////////////////////////////////

VtkDataWriter3D::VtkDataWriter3D(std::string const& fileName_, bool pointData_, bool mainProcOnly_)
    : fileName(fileName_),
      ostr(0),
      pointData(pointData_),
      mainProcOnly(mainProcOnly_)
{
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fileName.c_str());
        if (!(*ostr)) {
            std::cerr << "could not open file " <<  fileName << "\n";
            return;
        }
    }
}

VtkDataWriter3D::~VtkDataWriter3D() {
    delete ostr;
}

void VtkDataWriter3D::writeHeader(Box3D domain, Array<double,3> origin, double deltaX)
{
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        (*ostr) << "<?xml version=\"1.0\"?>\n";
#ifdef PLB_BIG_ENDIAN
        (*ostr) << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"BigEndian\">\n";
#else
        (*ostr) << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
#endif
        if (!pointData) {
            origin -= Array<double,3>(deltaX/2.,deltaX/2.,deltaX/2.);
            domain.x1++;
            domain.y1++;
            domain.z1++;
        }
        (*ostr) << "<ImageData WholeExtent=\""
                << domain.x0 << " " << domain.x1 << " "
                << domain.y0 << " " << domain.y1 << " "
                << domain.z0 << " " << domain.z1 << "\" "
                << "Origin=\""
                << origin[0] << " " << origin[1] << " " << origin[2] << "\" "
                << "Spacing=\""
                << deltaX << " " << deltaX << " " << deltaX << "\">\n";
    }
}

void VtkDataWriter3D::startPiece(Box3D domain) {
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        if (!pointData) {
            domain.x1++;
            domain.y1++;
            domain.z1++;
        }
        (*ostr) << "<Piece Extent=\""
                << domain.x0 << " " << domain.x1 << " "
                << domain.y0 << " " << domain.y1 << " "
                << domain.z0 << " " << domain.z1 << "\">\n";
        if (pointData) {
            (*ostr) << "<PointData>\n";
        }
        else {
            (*ostr) << "<CellData>\n";
        }
    }
}

void VtkDataWriter3D::endPiece() {
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        if (pointData) {
            (*ostr) << "</PointData>\n";
        }
        else {
            (*ostr) << "</CellData>\n";
        }
        (*ostr) << "</Piece>\n";
    }
}

void VtkDataWriter3D::writeFooter() {
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        (*ostr) << "</ImageData>\n";
        (*ostr) << "</VTKFile>\n";
    }
}


////////// class ParallelVtkDataWriter3D ////////////////////////////////////////

ParallelVtkDataWriter3D::ParallelVtkDataWriter3D(std::string const& fileName_)
    : sizeOfHeader(0),
      fileName(fileName_),
      ostr(0)
{ }

ParallelVtkDataWriter3D::~ParallelVtkDataWriter3D() { }

void ParallelVtkDataWriter3D::writeDataField( MultiBlock3D& block, IndexOrdering::OrderingT ordering, plint sizeOfScalar,
                                              plint dataSize)
{
    if (global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fileName.c_str(), std::ios_base::app|std::ios_base::binary);
        ostr->write((char const*)&dataSize, sizeof(dataSize));
        delete ostr;
    }
    global::mpi().barrier();
    parallelIO::saveFull(block, fileName, ordering, true);
    global::mpi().barrier();
}

void ParallelVtkDataWriter3D::writeHeader(Box3D domain, Array<double,3> origin, double deltaX)
{
    if (global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fileName.c_str(), std::ios_base::binary);
        if (!(*ostr)) {
            std::cerr << "could not open file " <<  fileName << "\n";
            return;
        }
        (*ostr) << "<?xml version=\"1.0\"?>\n";
#ifdef PLB_BIG_ENDIAN
        (*ostr) << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"BigEndian\" header_type=\"UInt64\">\n";
#else
        (*ostr) << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
#endif
        (*ostr) << "<ImageData WholeExtent=\""
                << domain.x0 << " " << domain.x1 << " "
                << domain.y0 << " " << domain.y1 << " "
                << domain.z0 << " " << domain.z1 << "\" "
                << "Origin=\""
                << origin[0] << " " << origin[1] << " " << origin[2] << "\" "
                << "Spacing=\""
                << deltaX << " " << deltaX << " " << deltaX << "\">\n";
        delete ostr;
    }
}

void ParallelVtkDataWriter3D::startPiece(Box3D domain) {
    if (global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fileName.c_str(), std::ios_base::app|std::ios_base::binary);
        if (!(*ostr)) {
            std::cerr << "could not open file " <<  fileName << "\n";
            return;
        }
        (*ostr) << "<Piece Extent=\""
                << domain.x0 << " " << domain.x1 << " "
                << domain.y0 << " " << domain.y1 << " "
                << domain.z0 << " " << domain.z1 << "\">\n";
        (*ostr) << "<PointData>\n";
        sizeOfHeader = ostr->tellp();
        delete ostr;
    }
}

void ParallelVtkDataWriter3D::endPiece() {
    if (global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fileName.c_str(), std::ios_base::app|std::ios_base::binary);
        if (!(*ostr)) {
            std::cerr << "could not open file " <<  fileName << "\n";
            return;
        }
        (*ostr) << "\n</PointData>\n";
        (*ostr) << "</Piece>\n";
        delete ostr;
    }
}

void ParallelVtkDataWriter3D::writeFooter1() {
    if (global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fileName.c_str(), std::ios_base::app|std::ios_base::binary);
        if (!(*ostr)) {
            std::cerr << "could not open file " <<  fileName << "\n";
        }
        (*ostr) << "</ImageData>\n";
        (*ostr) << "\n<AppendedData encoding=\"raw\">\n";
        (*ostr) << "   _";
        delete ostr;
    }
}

void ParallelVtkDataWriter3D::writeFooter2() {
    if (global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fileName.c_str(), std::ios_base::app|std::ios_base::binary);
        if (!(*ostr)) {
            std::cerr << "could not open file " <<  fileName << "\n";
            return;
        }
        (*ostr) << "\n</AppendedData>\n";
        (*ostr) << "\n</VTKFile>\n";
        delete ostr;
    }
}

void ParallelVtkDataWriter3D::appendFillerSpace(plint fillerSpace) {
    if (global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fileName.c_str(), std::ios_base::app|std::ios_base::binary);
        if (!(*ostr)) {
            std::cerr << "could not open file " <<  fileName << "\n";
            return;
        }
        std::string dummyString(fillerSpace, ' ');
        (*ostr) << dummyString;
        delete ostr;
    }
}




template<>
std::string VtkTypeNames<bool>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<char>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<unsigned char>::getBaseName() {
    return "UInt";
}

template<>
std::string VtkTypeNames<short int>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<unsigned short int>::getBaseName() {
    return "UInt";
}

template<>
std::string VtkTypeNames<int>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<unsigned int>::getBaseName() {
    return "UInt";
}

template<>
std::string VtkTypeNames<long int>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<unsigned long int>::getBaseName() {
    return "UInt";
}

template<>
std::string VtkTypeNames<float>::getBaseName() {
    return "Float";
}

template<>
std::string VtkTypeNames<double>::getBaseName() {
    return "Float";
}

template<>
std::string VtkTypeNames<long double>::getBaseName() {
    return "Float";
}

}  // namespace plb
