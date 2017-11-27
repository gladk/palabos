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
 * I/O routines for 2D multiblock -- header file.
 */
#ifndef MULTI_BLOCK_READER_2D_H
#define MULTI_BLOCK_READER_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/atomicContainerBlock2D.h"
#include "multiBlock/multiBlock2D.h"
#include "core/serializer.h"
#include "io/plbFiles.h"

namespace plb {

namespace parallelIO {

void dumpRestoreData( MultiBlock2D& multiBlock, bool dynamicContent,
                      std::vector<plint> const& myBlockIds, std::vector<std::vector<char> > const& data,
                      std::map<int,std::string> const& foreignIds = std::map<int,std::string>() );

void readkXmlSpec (
    FileName fName, Box2D& boundingBox, std::vector<plint>& offsets, plint& envelopeWidth,
    std::string& dataType, std::string& descriptor, std::string& family,
    std::vector<Box2D>& components, bool& dynamicContent, std::string& data_fName );

MultiBlock2D* load2D(FileName fName);

void load(FileName fName, MultiBlock2D& intoBlock, bool dynamicContent = true );

class SavedFullMultiBlockSerializer2D : public DataSerializer {
public:
    SavedFullMultiBlockSerializer2D(FileName fName);
    ~SavedFullMultiBlockSerializer2D();
    SavedFullMultiBlockSerializer2D(SavedFullMultiBlockSerializer2D const& rhs);
    virtual SavedFullMultiBlockSerializer2D* clone() const;
    virtual pluint getSize() const;
    virtual const char* getNextDataBuffer(pluint& bufferSize) const;
    virtual bool isEmpty() const;
    plint getCellDim() const;
    std::string dataType() const;
    Box2D getBoundingBox() const;
    bool orderingIsForward() const;
private:
    void computeSlice() const;
private:
    Box2D boundingBox;
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

#endif  // MULTI_BLOCK_READER_2D_H
