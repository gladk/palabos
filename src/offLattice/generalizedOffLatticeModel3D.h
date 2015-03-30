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

#ifndef GENERALIZED_OFF_LATTICE_MODEL_3D_H
#define GENERALIZED_OFF_LATTICE_MODEL_3D_H

#include "core/globalDefs.h"
#include "offLattice/offLatticeModel3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class ExtrapolatedGeneralizedOffLatticeModel3D : public OffLatticeModel3D<T,Array<T,3> >
{
public:
    ExtrapolatedGeneralizedOffLatticeModel3D(BoundaryShape3D<T,Array<T,3> >* shape_, int flowType_);
    ExtrapolatedGeneralizedOffLatticeModel3D(ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor> const& rhs);
    ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>& operator=(ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor> const& rhs);
    virtual ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>* clone() const;
    virtual plint getNumNeighbors() const;
    virtual void prepareCell (
            Dot3D const& cellLocation, AtomicContainerBlock3D& container );
    virtual void boundaryCompletion (
            AtomicBlock3D& lattice, AtomicContainerBlock3D& container,
            std::vector<AtomicBlock3D const*> const& args );
    virtual ContainerBlockData* generateOffLatticeInfo() const;
    virtual Array<T,3> getLocalForce(AtomicContainerBlock3D& container) const;
private:
    void cellCompletion (
            BlockLattice3D<T,Descriptor>& lattice,
            Dot3D const& genNode,
            std::vector<std::pair<int,int> > const& dryNodeFluidDirections,
            std::vector<int > const& dryNodeFluidNoSolidDirections,
            std::vector<plint> const& dryNodeIds, Dot3D const& absoluteOffset,
            Array<T,3>& localForce );
    void computeVelocity (
              BlockLattice3D<T,Descriptor> const& lattice, Dot3D const& genNode,
              Dot3D const& fluidDirection, int depth, Array<T,3> const& wallNode, T delta,
              Array<T,3> const& wall_vel, OffBoundary::Type bdType, Array<T,3> const& wallNormal,
              Array<T,Descriptor<T>::d>& u) const;
private:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class ExtrapolatedGeneralizedOffLatticeInfo3D : public ContainerBlockData {
    public:
        std::vector<Dot3D> const&                               getDryNodes() const
        { return dryNodes; }
        std::vector<Dot3D>&                                     getDryNodes()
        { return dryNodes; }
        std::vector<std::vector<std::pair<int,int> > > const&   getDryNodeFluidDirections() const
        { return dryNodeFluidDirections; }
        std::vector<std::vector<std::pair<int,int> > >&         getDryNodeFluidDirections()
        { return dryNodeFluidDirections; }
        std::vector<std::vector<int> > const&   getDryNodeFluidWithFluidDirections() const
        { return dryNodeFluidNoSolidDirections; }
        std::vector<std::vector<int> >&         getDryNodeFluidWithFluidDirections()
        { return dryNodeFluidNoSolidDirections; }
        std::vector<std::vector<plint> > const&                 getDryNodeIds() const
        { return dryNodeIds; }
        std::vector<std::vector<plint> >&                       getDryNodeIds()
        { return dryNodeIds; }
        Array<T,3> const&                                       getLocalForce() const
        { return localForce; }
        Array<T,3>&                                             getLocalForce()
        { return localForce; }
        virtual ExtrapolatedGeneralizedOffLatticeInfo3D* clone() const {
            return new ExtrapolatedGeneralizedOffLatticeInfo3D(*this);
        }
    private:
        std::vector<Dot3D> dryNodes;
        std::vector<std::vector<std::pair<int,int> > >   dryNodeFluidDirections;
        std::vector<std::vector<int> >   dryNodeFluidNoSolidDirections;
        std::vector<std::vector<plint> >                 dryNodeIds;
        Array<T,3>                                       localForce;
    };
};


template<typename T, template<typename U> class Descriptor>
class InterpolatedGeneralizedOffLatticeModel3D : public OffLatticeModel3D<T,Array<T,3> >
{
public:
    InterpolatedGeneralizedOffLatticeModel3D(BoundaryShape3D<T,Array<T,3> >* shape_, int flowType_);
    InterpolatedGeneralizedOffLatticeModel3D(InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor> const& rhs);
    InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>& operator=(InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor> const& rhs);
    virtual InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>* clone() const;
    virtual plint getNumNeighbors() const;
    virtual void prepareCell (
        Dot3D const& cellLocation, AtomicContainerBlock3D& container );
    virtual void boundaryCompletion (
        AtomicBlock3D& lattice, AtomicContainerBlock3D& container,
        std::vector<AtomicBlock3D const*> const& args );
    virtual ContainerBlockData* generateOffLatticeInfo() const;
    virtual Array<T,3> getLocalForce(AtomicContainerBlock3D& container) const;
private:
    void cellCompletion (
        BlockLattice3D<T,Descriptor>& lattice,
        Dot3D const& genNode,
        std::vector<std::pair<int,int> > const& wetNodeSolidDirections,
        std::vector<int > const& wetNodeFluidDirections,
        std::vector<plint> const& wetNodeIds, Dot3D const& absoluteOffset,
        Array<T,3>& localForce );
    void computeVelocity (
        BlockLattice3D<T,Descriptor> const& lattice, Dot3D const& genNode,
        Dot3D const& soldDirection, int depth, Array<T,3> const& wallNode, T wallDistance, T cellDistance,
        Array<T,3> const& wall_vel, OffBoundary::Type bdType, Array<T,3> const& wallNormal,
        Array<T,Descriptor<T>::d>& u) const;
private:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class InterpolatedGeneralizedOffLatticeInfo3D : public ContainerBlockData {
    public:
        std::vector<Dot3D> const&                               getWetNodes() const
        { return wetNodes; }
        std::vector<Dot3D>&                                     getWetNodes()
        { return wetNodes; }
        std::vector<std::vector<std::pair<int,int> > > const&   getWetNodeSolidDirections() const
        { return wetNodeSolidDirections; }
        std::vector<std::vector<std::pair<int,int> > >&         getWetNodeSolidDirections()
        { return wetNodeSolidDirections; }
        std::vector<std::vector<int> > const&   getWetNodeFluidDirections() const
        { return wetNodeFluidDirections; }
        std::vector<std::vector<int> >&         getWetNodeFluidDirections()
        { return wetNodeFluidDirections; }
        std::vector<std::vector<plint> > const&                 getWetNodeIds() const
        { return wetNodeIds; }
        std::vector<std::vector<plint> >&                       getWetNodeIds()
        { return wetNodeIds; }
        Array<T,3> const&                                       getLocalForce() const
        { return localForce; }
        Array<T,3>&                                             getLocalForce()
        { return localForce; }
        virtual InterpolatedGeneralizedOffLatticeInfo3D* clone() const {
            return new InterpolatedGeneralizedOffLatticeInfo3D(*this);
        }
    private:
        std::vector<Dot3D> wetNodes;
        std::vector<std::vector<std::pair<int,int> > > wetNodeSolidDirections; // stores the directions where there is a wall and valid neighbors in its opposite direction
        std::vector<std::vector<int> >                 wetNodeFluidDirections; // 
        std::vector<std::vector<plint> >               wetNodeIds;
        Array<T,3>                                     localForce;
    };
};

}  // namespace plb

#endif  // GUO_OFF_LATTICE_MODEL_3D_H
