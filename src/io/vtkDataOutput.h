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


#ifndef VTK_DATA_OUTPUT_H
#define VTK_DATA_OUTPUT_H

#include "core/globalDefs.h"
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "core/serializer.h"
#include "atomicBlock/dataField2D.h"
#include "multiBlock/multiDataField2D.h"
#include "atomicBlock/dataField3D.h"
#include "multiBlock/multiDataField3D.h"
#include "core/array.h"

namespace plb {

class VtkDataWriter3D {
public:
    VtkDataWriter3D(std::string const& fileName_);
    ~VtkDataWriter3D();
    void writeHeader(Box3D domain, Array<double,3> origin, double deltaX);
    void startPiece(Box3D domain);
    void endPiece();
    void writeFooter();
    template<typename T>
    void writeDataField( DataSerializer const* serializer,
                         std::string const& name, plint nDim );
private:
    VtkDataWriter3D(VtkDataWriter3D const& rhs);
    VtkDataWriter3D operator=(VtkDataWriter3D const& rhs);
private:
    std::string fileName;
    std::ofstream *ostr;
};

template<typename T>
class VtkImageOutput2D {
public:
    VtkImageOutput2D(std::string fName, double deltaX_=1.);
    VtkImageOutput2D(std::string fName, double deltaX_, Array<double,2> offset);
    ~VtkImageOutput2D();
    template<typename TConv>
    void writeData( plint nx, plint ny, plint nDim,
                    DataSerializer const* serializer, std::string const& name );
    template<typename TConv>
    void writeData(ScalarField2D<T>& scalarField,
                   std::string scalarFieldName, TConv scalingFactor=(T)1, TConv additiveOffset=(T)0);
    template<typename TConv>
    void writeData(MultiScalarField2D<T>& scalarField,
                   std::string scalarFieldName, TConv scalingFactor=(T)1, TConv additiveOffset=(T)0);
    template<plint n, typename TConv>
    void writeData(TensorField2D<T,n>& tensorField,
                   std::string tensorFieldName, TConv scalingFactor=(T)1);
    template<plint n, typename TConv>
    void writeData(MultiTensorField2D<T,n>& tensorField,
                   std::string tensorFieldName, TConv scalingFactor=(T)1);
    template<typename TConv>
    void writeData(MultiNTensorField2D<T>& nTensorField, std::string nTensorFieldName);
private:
    void writeHeader(plint nx_, plint ny_);
    void writeFooter();
private:
    std::string fullName;
    VtkDataWriter3D vtkOut;
    double deltaX;
    Array<T,2> offset;
    bool headerWritten;
    plint nx, ny;
};


template<typename T>
class VtkImageOutput3D {
public:
    VtkImageOutput3D(std::string fName, double deltaX_=1.);
    VtkImageOutput3D(std::string fName, double deltaX_, Array<double,3> offset);
    ~VtkImageOutput3D();
    template<typename TConv>
    void writeData( plint nx, plint ny, plint nz, plint nDim,
                    DataSerializer const* serializer, std::string const& name );
    template<typename TConv>
    void writeData( Box3D boundingBox_, plint nDim,
                    DataSerializer const* serializer, std::string const& name );
    template<typename TConv>
    void writeData(ScalarField3D<T>& scalarField,
                   std::string scalarFieldName, TConv scalingFactor=(T)1, TConv additiveOffset=(T)0);
    template<typename TConv>
    void writeData(MultiScalarField3D<T>& scalarField,
                   std::string scalarFieldName, TConv scalingFactor=(T)1, TConv additiveOffset=(T)0);
    template<plint n, typename TConv>
    void writeData(TensorField3D<T,n>& tensorField,
                   std::string tensorFieldName, TConv scalingFactor=(T)1);
    template<plint n, typename TConv>
    void writeData(MultiTensorField3D<T,n>& tensorField,
                   std::string tensorFieldName, TConv scalingFactor=(T)1);
    template<typename TConv>
    void writeData(MultiNTensorField3D<T>& nTensorField, std::string nTensorFieldName);
private:
    void writeHeader(plint nx_, plint ny_, plint nz_);
    void writeHeader(Box3D boundingBox_);
    void writeFooter();
private:
    std::string fullName;
    VtkDataWriter3D vtkOut;
    double deltaX;
    Array<T,3> offset;
    bool headerWritten;
    Box3D boundingBox;
};

} // namespace plb

#endif  // VTK_DATA_OUTPUT_H
