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
 * Serial implementation of scalar, vector and tensor fields for 3D data analysis.
 * -- header file
 */

#ifndef MULTI_CONTAINER_BLOCK_3D_H
#define MULTI_CONTAINER_BLOCK_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/atomicContainerBlock3D.h"
#include "multiBlock/multiBlock3D.h"

namespace plb {

class MultiContainerBlock3D : public MultiBlock3D {
public:
    typedef std::map<plint,AtomicContainerBlock3D*> BlockMap;
public:
    MultiContainerBlock3D (
            MultiBlockManagement3D const& multiBlockManagement_,
            CombinedStatistics* combinedStatistics_ );
    MultiContainerBlock3D(plint nx_, plint ny_, plint nz_);
    MultiContainerBlock3D(MultiBlock3D const& rhs);
    MultiContainerBlock3D(MultiBlock3D const& rhs, Box3D subDomain, bool crop);
    ~MultiContainerBlock3D();
    MultiContainerBlock3D& operator=(MultiContainerBlock3D const& rhs);
    MultiContainerBlock3D(MultiContainerBlock3D const& rhs);
    MultiContainerBlock3D* clone() const;
    MultiContainerBlock3D* clone(MultiBlockManagement3D const& multiBlockManagement) const;
    void swap(MultiContainerBlock3D& rhs);
public:
    virtual AtomicContainerBlock3D& getComponent(plint iBlock);
    virtual AtomicContainerBlock3D const& getComponent(plint iBlock) const;
    virtual plint sizeOfCell() const;
    virtual plint getCellDim() const;
    virtual int getStaticId() const;
    virtual void copyReceive (
                MultiBlock3D const& fromBlock, Box3D const& fromDomain,
                Box3D const& toDomain, modif::ModifT whichData=modif::dataStructure );
    std::string getBlockName() const;
    std::vector<std::string> getTypeInfo() const;
private:
    void allocateBlocks();
    void allocateBlocks(MultiContainerBlock3D const& rhs);
    void deAllocateBlocks();
private:
    BlockMap blocks;
};

MultiContainerBlock3D* createContainerBlock(MultiBlock3D& templ, ContainerBlockData* data);

}  // namespace plb

#endif  // MULTI_CONTAINER_BLOCK_3D_H
