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
 * Helper classes for parallel 3D multiblock lattice -- header file.
 */

#ifndef PARALLEL_BLOCK_COMMUNICATOR_3D_H
#define PARALLEL_BLOCK_COMMUNICATOR_3D_H

#include "core/globalDefs.h"
#include "multiBlock/blockCommunicator3D.h"
#include "multiBlock/multiBlockManagement3D.h"
#include "multiBlock/multiBlock3D.h"
#include "parallelism/sendRecvPool.h"
#include "parallelism/communicationPackage3D.h"
#include <vector>

namespace plb {

#ifdef PLB_MPI_PARALLEL

struct CommunicationPattern3D
{
    CommunicationPattern3D (
            std::vector<Overlap3D> const& overlaps,
            MultiBlockManagement3D const& originManagement,
            MultiBlockManagement3D const& destinationManagement );
    CommunicationPackage3D sendPackage;
    CommunicationPackage3D recvPackage;
    CommunicationPackage3D sendRecvPackage;
};

struct CommunicationStructure3D
{
    CommunicationStructure3D (
            std::vector<Overlap3D> const& overlaps,
            MultiBlockManagement3D const& originManagement,
            MultiBlockManagement3D const& destinationManagement,
            plint sizeOfCell );
    CommunicationPackage3D sendPackage;
    CommunicationPackage3D recvPackage;
    CommunicationPackage3D sendRecvPackage;
    SendPoolCommunicator sendComm;
    RecvPoolCommunicator recvComm;
};


class ParallelBlockCommunicator3D : public BlockCommunicator3D {
public:
    ParallelBlockCommunicator3D();
    ~ParallelBlockCommunicator3D();
    ParallelBlockCommunicator3D(ParallelBlockCommunicator3D const& rhs);
    ParallelBlockCommunicator3D& operator= (
            ParallelBlockCommunicator3D const& rhs );
    void swap(ParallelBlockCommunicator3D& rhs);
    virtual ParallelBlockCommunicator3D* clone() const;
    virtual void duplicateOverlaps(MultiBlock3D& multiBlock, modif::ModifT whichData) const;
    virtual void communicate( std::vector<Overlap3D> const& overlaps,
                              MultiBlock3D const& originMultiBlock,
                              MultiBlock3D& destinationMultiBlock,
                              modif::ModifT whichData ) const;
    virtual void signalPeriodicity() const;
private:
    void communicate( CommunicationStructure3D& communication,
                      MultiBlock3D const& originMultiBlock,
                      MultiBlock3D& destinationMultiBlock, modif::ModifT whichData ) const;
    void subscribeOverlap (
        Overlap3D const& overlap, MultiBlockManagement3D const& multiBlockManagement,
        SendRecvPool& sendPool, SendRecvPool& recvPool, plint sizeOfCell ) const;
private:
    mutable bool overlapsModified;
    mutable CommunicationStructure3D* communication;
};


class BlockingCommunicator3D : public BlockCommunicator3D {
public:
    BlockingCommunicator3D();
    ~BlockingCommunicator3D();
    BlockingCommunicator3D(BlockingCommunicator3D const& rhs);
    BlockingCommunicator3D& operator= (
            BlockingCommunicator3D const& rhs );
    void swap(BlockingCommunicator3D& rhs);
    virtual BlockingCommunicator3D* clone() const;
    virtual void duplicateOverlaps(MultiBlock3D& multiBlock, modif::ModifT whichData) const;
    virtual void communicate( std::vector<Overlap3D> const& overlaps,
                              MultiBlock3D const& originMultiBlock,
                              MultiBlock3D& destinationMultiBlock,
                              modif::ModifT whichData ) const;
    virtual void signalPeriodicity() const;
private:
    void communicate( CommunicationPattern3D& communication,
                      MultiBlock3D const& originMultiBlock,
                      MultiBlock3D& destinationMultiBlock, modif::ModifT whichData ) const;
private:
    mutable bool overlapsModified;
    mutable CommunicationPattern3D* communication;
};

#endif  // PLB_MPI_PARALLEL

}  // namespace plb

#endif  // PARALLEL_BLOCK_COMMUNICATOR_3D_H
