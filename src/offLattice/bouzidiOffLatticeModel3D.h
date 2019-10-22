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

#ifndef BOUZIDI_OFF_LATTICE_MODEL_3D_H
#define BOUZIDI_OFF_LATTICE_MODEL_3D_H

#include "core/globalDefs.h"
#include "offLattice/offLatticeModel3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class BouzidiOffLatticeModel3D : public OffLatticeModel3D<T,Array<T,3> >
{
public:
    BouzidiOffLatticeModel3D(BoundaryShape3D<T,Array<T,3> >* shape_, int flowType_);
    virtual BouzidiOffLatticeModel3D<T,Descriptor>* clone() const;
    virtual plint getNumNeighbors() const;
    virtual bool isExtrapolated() const;
    virtual void prepareCell (
            Dot3D const& cellLocation, AtomicContainerBlock3D& container );
    virtual void boundaryCompletion (
            AtomicBlock3D& lattice, AtomicContainerBlock3D& container,
            std::vector<AtomicBlock3D *> const& args );
    void cellCompletion (
            BlockLattice3D<T,Descriptor>& lattice,
            Dot3D const& boundaryNode,
            std::vector<int> const& solidDirections, std::vector<plint> const& boundaryIds,
            std::vector<bool> const& hasFluidNeighbor, Dot3D const& absoluteOffset,
            Array<T,3>& localForce, std::vector<AtomicBlock3D *> const& args );
    virtual ContainerBlockData* generateOffLatticeInfo() const;
    virtual Array<T,3> getLocalForce(AtomicContainerBlock3D& container) const;
private:
    std::vector<T> invAB;
private:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class BouzidiOffLatticeInfo3D : public ContainerBlockData {
    public:
        std::vector<Dot3D> const&               getBoundaryNodes() const
        { return boundaryNodes; }
        std::vector<Dot3D>&                     getBoundaryNodes()
        { return boundaryNodes; }
        std::vector<std::vector<int> > const&   getSolidDirections() const
        { return solidDirections; }
        std::vector<std::vector<int> >&         getSolidDirections()
        { return solidDirections; }
        std::vector<std::vector<plint> > const& getBoundaryIds() const
        { return boundaryIds; }
        std::vector<std::vector<plint> >&       getBoundaryIds()
        { return boundaryIds; }
        std::vector<std::vector<bool> > const&  getHasFluidNeighbor() const
        { return hasFluidNeighbor; }
        std::vector<std::vector<bool> >&        getHasFluidNeighbor()
        { return hasFluidNeighbor; }
        Array<T,3> const&                       getLocalForce() const
        { return localForce; }
        Array<T,3>&                             getLocalForce()
        { return localForce; }
        virtual BouzidiOffLatticeInfo3D* clone() const {
            return new BouzidiOffLatticeInfo3D(*this);
        }
    private:
        std::vector<Dot3D>               boundaryNodes;
        std::vector<std::vector<int> >   solidDirections;
        std::vector<std::vector<plint> > boundaryIds;
        std::vector<std::vector<bool> >  hasFluidNeighbor;
        Array<T,3>                       localForce;
    };
};

}  // namespace plb

#endif  // BOUZIDI_OFF_LATTICE_MODEL_3D_H

