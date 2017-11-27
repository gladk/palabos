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

/** \file
 * I/O routines for 3D multiblock -- header file.
 */
#ifndef MULTI_BLOCK_READER_3D_H
#define MULTI_BLOCK_READER_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/atomicContainerBlock3D.h"
#include "multiBlock/multiBlock3D.h"
#include "core/serializer.h"
#include "io/plbFiles.h"

namespace plb {

namespace parallelIO {

void dumpRestoreData( MultiBlock3D& multiBlock, bool dynamicContent,
                      std::vector<plint> const& myBlockIds, std::vector<std::vector<char> > const& data,
                      std::map<int,std::string> const& foreignIds = std::map<int,std::string>() );

void readkXmlSpec (
    FileName fName, Box3D& boundingBox, std::vector<plint>& offsets, plint& envelopeWidth,
    std::string& dataType, std::string& descriptor, std::string& family,
    std::vector<Box3D>& components, bool& dynamicContent, std::string& data_fName );

MultiBlock3D* load3D(FileName fName);

void load(FileName fName, MultiBlock3D& intoBlock, bool dynamicContent = true );

class SavedFullMultiBlockSerializer3D : public DataSerializer {
public:
    SavedFullMultiBlockSerializer3D(FileName fName);
    ~SavedFullMultiBlockSerializer3D();
    SavedFullMultiBlockSerializer3D(SavedFullMultiBlockSerializer3D const& rhs);
    virtual SavedFullMultiBlockSerializer3D* clone() const;
    virtual pluint getSize() const;
    virtual const char* getNextDataBuffer(pluint& bufferSize) const;
    virtual bool isEmpty() const;
    plint getCellDim() const;
    std::string dataType() const;
    Box3D getBoundingBox() const;
    bool orderingIsForward() const;
private:
    void computeSlice() const;
private:
    Box3D boundingBox;
    plint cellDim, typeSize;
    plint sizeOfChunk;
    FileName data_fName;
    std::string str_dataType;
    bool forwardOrdering;
    mutable plint pos;
    mutable std::vector<char> buffer;
    FILE *fp;
};

}  // namespace parallelIO

}  // namespace plb

#endif  // MULTI_BLOCK_READER_3D_H
