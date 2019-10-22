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
 * Serializer and UnSerializer for multi blocks -- header file.
 */
#ifndef MULTI_BLOCK_SERIALIZER_2D_H
#define MULTI_BLOCK_SERIALIZER_2D_H

#include "core/globalDefs.h"
#include "multiBlock/multiBlock2D.h"
#include "core/serializer.h"
#include "core/util.h"

namespace plb {

class MultiBlockSerializer2D : public DataSerializer {
public:
    MultiBlockSerializer2D(MultiBlock2D const& multiBlock_,
                           IndexOrdering::OrderingT ordering_);
    MultiBlockSerializer2D(MultiBlock2D const& multiBlock_,
                           Box2D domain_,
                           IndexOrdering::OrderingT ordering_);
    virtual MultiBlockSerializer2D* clone() const;
    virtual pluint getSize() const;
    virtual const char* getNextDataBuffer(pluint& bufferSize) const;
    virtual bool isEmpty() const;
private:
    SparseBlockStructure2D const& getSparseBlockStructure() const;
    bool isLocal(plint blockId) const;
    void computeBufferAlongX(plint nextBlockId, plint nextChunkSize) const;
    void computeBufferAlongY(plint nextBlockId, plint nextChunkSize) const;
    void communicateBuffer(plint bufferSize, plint fromBlockId, bool isAllocated) const;
    void fillBufferWithZeros(plint nextChunkSize) const;
private:
    MultiBlock2D const& multiBlock;
    IndexOrdering::OrderingT ordering;
    Box2D domain;
    mutable plint iX, iY;
    mutable std::vector<char> buffer;
};

class MultiBlockUnSerializer2D : public DataUnSerializer {
public:
    MultiBlockUnSerializer2D(MultiBlock2D& multiBlock_,
                             IndexOrdering::OrderingT ordering_);
    MultiBlockUnSerializer2D(MultiBlock2D& multiBlock_,
                             Box2D domain_,
                             IndexOrdering::OrderingT ordering_);
    virtual MultiBlockUnSerializer2D* clone() const;
    virtual pluint getSize() const;
    virtual char* getNextDataBuffer(pluint& bufferSize);
    virtual void commitData();
    virtual bool isFull() const;
private:
    SparseBlockStructure2D const& getSparseBlockStructure() const;
    bool isLocal(plint blockId) const;
    void fillBufferAlongX(plint nextBlockId, plint nextChunkSize);
    void fillBufferAlongY(plint nextBlockId, plint nextChunkSize);
    void communicateBuffer(plint bufferSize, plint toBlockId, bool isAllocated) const;
private:
    MultiBlock2D& multiBlock;
    IndexOrdering::OrderingT ordering;
    Box2D domain;
    mutable plint iX, iY;
    mutable std::vector<char> buffer;
};

class MultiBlockFastSerializer2D : public DataSerializer {
public:
    MultiBlockFastSerializer2D(MultiBlock2D const& multiBlock_,
                               IndexOrdering::OrderingT ordering_);
    MultiBlockFastSerializer2D(MultiBlock2D const& multiBlock_,
                               Box2D domain_,
                               IndexOrdering::OrderingT ordering_);
    virtual MultiBlockFastSerializer2D* clone() const;
    virtual pluint getSize() const;
    virtual const char* getNextDataBuffer(pluint& bufferSize) const;
    virtual bool isEmpty() const;
private:
    pluint computeSlice() const;
private:
    MultiBlock2D const& multiBlock;
    IndexOrdering::OrderingT ordering;
    Box2D domain;
    mutable plint pos;
    mutable std::vector<char> buffer;
};

class MultiBlockFastUnSerializer2D : public DataUnSerializer {
public:
    MultiBlockFastUnSerializer2D(MultiBlock2D& multiBlock_,
                                 IndexOrdering::OrderingT ordering_);
    MultiBlockFastUnSerializer2D(MultiBlock2D& multiBlock_,
                                 Box2D domain_,
                                 IndexOrdering::OrderingT ordering_);
    virtual MultiBlockFastUnSerializer2D* clone() const;
    virtual pluint getSize() const;
    virtual char* getNextDataBuffer(pluint& bufferSize);
    virtual void commitData();
    virtual bool isFull() const;
private:
    MultiBlock2D& multiBlock;
    IndexOrdering::OrderingT ordering;
    Box2D domain;
    mutable plint pos;
    mutable std::vector<char> buffer;
};
}  //  namespace plb

#endif  // MULTI_BLOCK_SERIALIZER_2D_H
