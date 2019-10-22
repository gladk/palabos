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


#ifndef SPARSE_VTK_DATA_OUTPUT_H
#define SPARSE_VTK_DATA_OUTPUT_H

#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "core/globalDefs.h"
#include "io/plbFiles.h"
#include "io/vtkDataOutput.h"
#include "multiBlock/group3D.h"

namespace plb {

class SparseVtkImageOutput3D {
public:
    SparseVtkImageOutput3D(FileName fName_);
    ~SparseVtkImageOutput3D();
    void writeVtkBlock( Group3D& group, double deltaX, Array<double,3> const& offset,
                        plint vtkBlockId=0, bool pointData=true );
    void writeVtkBlock( Group3D& group, double deltaX, plint vtkBlockId=0, bool pointData=true );
private:
    std::string writeAtomicBlock (
            Group3D& group, plint partId, plint vtkBlockId, plint atomicBlockId,
            Box3D bulk, bool pointData, bool isLocal, double deltaX, Array<double,3> offset );
    void writeField (
            MultiBlock3D const& multiBlock, plint atomicBlockId, Box3D bulk,
            std::string fieldName, VtkDataWriter3D& vtkOut );
private:
    FileName fName;
    std::string dName, rdName;
    plb_ofstream vtmFile;
    Array<double,3> offset;
};


} // namespace plb

#endif  // SPARSE_VTK_DATA_OUTPUT_H
