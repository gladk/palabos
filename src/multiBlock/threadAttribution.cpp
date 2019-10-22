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
 * Attribution of block id to thread for parallel execution -- implementation.
 */

#include "parallelism/mpiManager.h"
#include "multiBlock/threadAttribution.h"
#include "core/plbDebug.h"

namespace plb {

bool SerialThreadAttribution::isLocal(plint blockId) const {
    return true;
}

bool SerialThreadAttribution::allBlocksAreLocal() const {
    return true;
}

int SerialThreadAttribution::getMpiProcess(plint blockId) const {
    return global::mpi().bossId();
}

int SerialThreadAttribution::getLocalThreadId(plint blockId) const {
    return 0;
}

ThreadAttribution* SerialThreadAttribution::merge (
            ThreadAttribution const& rhs,
            std::map<plint,std::vector<plint> > const& remappedIds ) const
{
    // Serial remains serial. Nothing to be done here.
    return new SerialThreadAttribution();
}

ThreadAttribution* SerialThreadAttribution::extend (
              std::vector<plint> const& ids,
              std::vector<plint> const& mpiProcesses, std::vector<plint> const& localThreads ) const
{
    // Serial remains serial. Nothing to be done here.
    return new SerialThreadAttribution();
}

ThreadAttribution* SerialThreadAttribution::clone() const {
    return new SerialThreadAttribution;
}


bool OneToOneThreadAttribution::isLocal(plint blockId) const {
    return blockId == global::mpi().getRank();
}

bool OneToOneThreadAttribution::allBlocksAreLocal() const {
    return false;
}

int OneToOneThreadAttribution::getMpiProcess(plint blockId) const {
    return blockId;
}
int OneToOneThreadAttribution::getLocalThreadId(plint blockId) const {
    return 0;
}


ThreadAttribution* OneToOneThreadAttribution::merge (
            ThreadAttribution const& rhs,
            std::map<plint,std::vector<plint> > const& remappedIds ) const
{
    // Need to turn to an explicit form of thread attribution. Remember that
    //   ExplicitThreadAttribution defaults to OneToOneThreadAttribution
    //   for all non-specified blockIds, so we are safe.
    ExplicitThreadAttribution* newThreadAttribution
        = new ExplicitThreadAttribution();
    std::map<plint,std::vector<plint> >::const_iterator it = remappedIds.begin();
    for (; it != remappedIds.end(); ++it) {
        plint oldId = it->first;
        std::vector<plint> const& newIds = it->second;
        for (pluint iNew=0; iNew<newIds.size(); ++iNew) {
            newThreadAttribution->addBlock (
                    newIds[iNew],
                    rhs.getMpiProcess(oldId),
                    rhs.getLocalThreadId(oldId) );
        }
    }
    return newThreadAttribution;
}

ThreadAttribution* OneToOneThreadAttribution::extend (
              std::vector<plint> const& ids,
              std::vector<plint> const& mpiProcesses, std::vector<plint> const& localThreads ) const
{
    PLB_PRECONDITION( ids.size() == mpiProcesses.size() &&
                      ids.size() == localThreads.size() );
    // Need to turn to an explicit form of thread attribution. Remember that
    //   ExplicitThreadAttribution defaults to OneToOneThreadAttribution
    //   for all non-specified blockIds, so we are safe.
    ExplicitThreadAttribution* newThreadAttribution
        = new ExplicitThreadAttribution();
    for (pluint iNew=0; iNew<ids.size(); ++iNew) {
        newThreadAttribution->addBlock (
                ids[iNew], mpiProcesses[iNew], localThreads[iNew] );
    }
    return newThreadAttribution;
}

ThreadAttribution* OneToOneThreadAttribution::clone() const {
    return new OneToOneThreadAttribution;
}

ExplicitThreadAttribution::ExplicitThreadAttribution()
{ }

ExplicitThreadAttribution::ExplicitThreadAttribution (
        std::map<plint,plint> const& mpiProcessAttribution_ )
    : mpiProcessAttribution(mpiProcessAttribution_)
{
    std::map<plint,plint>::const_iterator it = mpiProcessAttribution.begin();
    // Default local threads to zero.
    for (; it != mpiProcessAttribution.end(); ++it) {
        localThreadAttribution[it->first] = 0;
    }
}

ExplicitThreadAttribution::ExplicitThreadAttribution (
        std::map<plint,plint> const& mpiProcessAttribution_,
        std::map<plint,plint> const& localThreadAttribution_ )
    : mpiProcessAttribution(mpiProcessAttribution_),
      localThreadAttribution(localThreadAttribution_)
{ }

void ExplicitThreadAttribution::addBlock (
        plint blockId, plint mpiProcess, plint localThread )
{
    mpiProcessAttribution[blockId] = mpiProcess;
    localThreadAttribution[blockId] = localThread;
}

bool ExplicitThreadAttribution::isLocal(plint blockId) const {
    std::map<plint,plint>::const_iterator it = mpiProcessAttribution.find(blockId);
    // If this blockId is not registered, fall back to default behavior,
    //   which is the one specified by OneToOneThreadAttribution.
    if (it == mpiProcessAttribution.end()) {
        return blockId == global::mpi().getRank();
    }
    else {
        return it->second == global::mpi().getRank();
    }
}

bool ExplicitThreadAttribution::allBlocksAreLocal() const {
    // They might all be local, but we a priori don't know. Remember that
    //   ExplicitThreadAttribution defaults to OneToOneThreadAttribution
    //   for all blockIds which have not been subscribed. Does, we can't
    //   know if all blocks are local or not.
    return false;
}

int ExplicitThreadAttribution::getMpiProcess(plint blockId) const {
    std::map<plint,plint>::const_iterator it = mpiProcessAttribution.find(blockId);
    // If this blockId is not registered, fall back to default behavior,
    //   which is the one specified by OneToOneThreadAttribution.
    if (it == mpiProcessAttribution.end()) {
        return blockId;
    }
    else {
        return it->second;
    }
}

int ExplicitThreadAttribution::getLocalThreadId(plint blockId) const {
    std::map<plint,plint>::const_iterator it = localThreadAttribution.find(blockId);
    // If this blockId is not registered, fall back to default behavior,
    //   which is the one specified by OneToOneThreadAttribution.
    if (it == localThreadAttribution.end()) {
        return 0;
    }
    else {
        return it->second;
    }
}

ThreadAttribution* ExplicitThreadAttribution::merge (
            ThreadAttribution const& rhs,
            std::map<plint,std::vector<plint> > const& remappedIds ) const
{
    ExplicitThreadAttribution* newThreadAttribution
        = new ExplicitThreadAttribution();
    newThreadAttribution->mpiProcessAttribution = mpiProcessAttribution;
    newThreadAttribution->localThreadAttribution = localThreadAttribution;

    std::map<plint,std::vector<plint> >::const_iterator it = remappedIds.begin();
    for (; it != remappedIds.end(); ++it) {
        plint oldId = it->first;
        std::vector<plint> const& newIds = it->second;
        for (pluint iNew=0; iNew<newIds.size(); ++iNew) {
            newThreadAttribution->addBlock (
                    newIds[iNew],
                    rhs.getMpiProcess(oldId),
                    rhs.getLocalThreadId(oldId) );
        }
    }
    return newThreadAttribution;
}

ThreadAttribution* ExplicitThreadAttribution::extend (
              std::vector<plint> const& ids,
              std::vector<plint> const& mpiProcesses, std::vector<plint> const& localThreads ) const
{
    PLB_PRECONDITION( ids.size() == mpiProcesses.size() &&
                      ids.size() == localThreads.size() );
    ExplicitThreadAttribution* newThreadAttribution
        = new ExplicitThreadAttribution(mpiProcessAttribution, localThreadAttribution);
    for (pluint iNew=0; iNew<ids.size(); ++iNew) {
        newThreadAttribution->addBlock (
                ids[iNew], mpiProcesses[iNew], localThreads[iNew] );
    }
    return newThreadAttribution;
}

ThreadAttribution* ExplicitThreadAttribution::clone() const {
    return new ExplicitThreadAttribution(*this);
}

}  // namespace plb
