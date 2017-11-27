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

#ifndef FILIPPOVA_HAENEL_OFF_LATTICE_MODEL_3D_H
#define FILIPPOVA_HAENEL_OFF_LATTICE_MODEL_3D_H

#include "core/globalDefs.h"
#include "offLattice/offLatticeModel3D.h"
#include "offLattice/guoOffLatticeModel3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class FilippovaHaenelModel3D : public OffLatticeModel3D<T,Array<T,3> >
{
public:
    FilippovaHaenelModel3D(BoundaryShape3D<T,Array<T,3> >* shape_, int flowType_, bool useAllDirections_=true);
    virtual FilippovaHaenelModel3D<T,Descriptor>* clone() const;
    virtual plint getNumNeighbors() const;
    virtual bool isExtrapolated() const;
    virtual void prepareCell (
            Dot3D const& cellLocation, AtomicContainerBlock3D& container );
    virtual void boundaryCompletion (
            AtomicBlock3D& lattice, AtomicContainerBlock3D& container,
            std::vector<AtomicBlock3D *> const& args );

    virtual ContainerBlockData* generateOffLatticeInfo() const;
    virtual Array<T,3> getLocalForce(AtomicContainerBlock3D& container) const;
    void selectComputeStat(bool flag) { computeStat = flag; }
    bool computesStat() const { return computeStat; }
private:
    void cellCompletion (
            BlockLattice3D<T,Descriptor>& lattice,
            Dot3D const& guoNode,
            std::vector<int> const& dryNodeFluidDirections,
            std::vector<plint> const& dryNodeIds, Dot3D const& absoluteOffset,
            Array<T,3>& localForce, std::vector<AtomicBlock3D *> const& args );
private:
    bool computeStat;
private:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class OffLatticeInfo3D : public ContainerBlockData {
    public:
        std::vector<Dot3D> const&                getDryNodes() const
        { return dryNodes; }
        std::vector<Dot3D>&                      getDryNodes()
        { return dryNodes; }
        std::vector<std::vector<int > > const&   getDryNodeFluidDirections() const
        { return dryNodeFluidDirections; }
        std::vector<std::vector<int > >&         getDryNodeFluidDirections()
        { return dryNodeFluidDirections; }
        std::vector<std::vector<plint> > const&  getDryNodeIds() const
        { return dryNodeIds; }
        std::vector<std::vector<plint> >&        getDryNodeIds()
        { return dryNodeIds; }
        Array<T,3> const&                        getLocalForce() const
        { return localForce; }
        Array<T,3>&                              getLocalForce()
        { return localForce; }
        virtual OffLatticeInfo3D* clone() const {
            return new OffLatticeInfo3D(*this);
        }
    private:
        std::vector<Dot3D> dryNodes;
        std::vector<std::vector<int> >   dryNodeFluidDirections;
        std::vector<std::vector<plint> > dryNodeIds;
        Array<T,3>                       localForce;
    };
};

}  // namespace plb

#endif  // FILIPPOVA_HAENEL_OFF_LATTICE_MODEL_3D_H

