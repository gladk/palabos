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

#ifndef GUO_OFF_LATTICE_MODEL_3D_H
#define GUO_OFF_LATTICE_MODEL_3D_H

#include "core/globalDefs.h"
#include "offLattice/offLatticeModel3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class GuoOffLatticeModel3D : public OffLatticeModel3D<T,Array<T,3> >
{
public:
    GuoOffLatticeModel3D(BoundaryShape3D<T,Array<T,3> >* shape_, int flowType_, bool useAllDirections_=true);
    virtual GuoOffLatticeModel3D<T,Descriptor>* clone() const;
    virtual plint getNumNeighbors() const;
    virtual bool isExtrapolated() const;
    virtual void prepareCell (
            Dot3D const& cellLocation, AtomicContainerBlock3D& container );
    virtual void boundaryCompletion (
            AtomicBlock3D& lattice, AtomicContainerBlock3D& container,
            std::vector<AtomicBlock3D *> const& args );
    virtual ContainerBlockData* generateOffLatticeInfo() const;
    virtual Array<T,3> getLocalForce(AtomicContainerBlock3D& container) const;
private:
    void cellCompletion (
            BlockLattice3D<T,Descriptor>& lattice,
            Dot3D const& guoNode,
            std::vector<std::pair<int,int> > const& dryNodeFluidDirections,
            std::vector<plint> const& dryNodeIds, Dot3D const& absoluteOffset,
            Array<T,3>& localForce, std::vector<AtomicBlock3D *> const& args );
    void computeRhoBarJPiNeqAlongDirection (
              BlockLattice3D<T,Descriptor> const& lattice, Dot3D const& guoNode,
              Dot3D const& fluidDirection, int depth, Array<T,3> const& wallNode, T delta,
              Array<T,3> const& wall_vel, OffBoundary::Type bdType,
              Array<T,3> const& wallNormal, plint triangleId,
              T& rhoBar, Array<T,Descriptor<T>::d>& j,
              Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq,
              std::vector<AtomicBlock3D *> const& args ) const;
    void computeRhoBarJfNeqAlongDirection (
              BlockLattice3D<T,Descriptor> const& lattice, Dot3D const& guoNode,
              Dot3D const& fluidDirection, int depth, Array<T,3> const& wallNode, T delta,
              Array<T,3> const& wall_vel, OffBoundary::Type bdType,
              Array<T,3> const& wallNormal, plint triangleId,
              T& rhoBar, Array<T,Descriptor<T>::d>& j,
              Array<T,Descriptor<T>::q>& fNeq,
              std::vector<AtomicBlock3D *> const& args ) const;
private:
    bool useAllDirections;
public:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class GuoOffLatticeInfo3D : public ContainerBlockData {
    public:
        std::vector<Dot3D> const&                               getDryNodes() const
        { return dryNodes; }
        std::vector<Dot3D>&                                     getDryNodes()
        { return dryNodes; }
        std::vector<std::vector<std::pair<int,int> > > const&   getDryNodeFluidDirections() const
        { return dryNodeFluidDirections; }
        std::vector<std::vector<std::pair<int,int> > >&         getDryNodeFluidDirections()
        { return dryNodeFluidDirections; }
        std::vector<std::vector<plint> > const&                 getDryNodeIds() const
        { return dryNodeIds; }
        std::vector<std::vector<plint> >&                       getDryNodeIds()
        { return dryNodeIds; }
        std::vector<bool> const&                                getIsConnected() const
        { return isConnected; }
        std::vector<bool>&                                      getIsConnected()
        { return isConnected; }
        Array<T,3> const&                                       getLocalForce() const
        { return localForce; }
        Array<T,3>&                                             getLocalForce()
        { return localForce; }
        virtual GuoOffLatticeInfo3D* clone() const {
            return new GuoOffLatticeInfo3D(*this);
        }
    private:
        std::vector<Dot3D> dryNodes;
        std::vector<std::vector<std::pair<int,int> > >   dryNodeFluidDirections;
        std::vector<std::vector<plint> >                 dryNodeIds;
        std::vector<bool>                                isConnected;
        Array<T,3>                                       localForce;
    };

    struct LiquidNeighbor {
        LiquidNeighbor(plint iNeighbor_, plint depth_, plint iTriangle_, Array<T,3> wallNormal);
        bool operator<(LiquidNeighbor const& rhs) const;
        plint iNeighbor, depth;
        plint iTriangle;
        T cosAngle;
    };
};

template<typename T, template<typename U> class Descriptor>
class GuoOffLatticeFdModel3D : public OffLatticeModel3D<T,Array<T,3> >
{
public:
    struct LiquidNeighbor {
        LiquidNeighbor(plint iNeighbor_, plint depth_, plint iTriangle_, Array<T,3> wallNormal);
        bool operator<(LiquidNeighbor const& rhs) const;
        plint iNeighbor, depth;
        plint iTriangle;
        T cosAngle;
    };
public:
    GuoOffLatticeFdModel3D(BoundaryShape3D<T,Array<T,3> >* shape_, int flowType_, bool useAllDirections_=true);
    virtual GuoOffLatticeFdModel3D<T,Descriptor>* clone() const;
    virtual plint getNumNeighbors() const;
    virtual bool isExtrapolated() const;
    virtual void prepareCell (
            Dot3D const& cellLocation, AtomicContainerBlock3D& container );
    virtual void boundaryCompletion (
            AtomicBlock3D& lattice, AtomicContainerBlock3D& container,
            std::vector<AtomicBlock3D *> const& args );
    virtual ContainerBlockData* generateOffLatticeInfo() const;
    virtual Array<T,3> getLocalForce(AtomicContainerBlock3D& container) const;
private:
    bool isUsable(const Dot3D &pos) const;
    std::pair<int,int> computeOrderAndDirection(const Dot3D &pos, const Dot3D &dx) const;
    void cellCompletion (
            BlockLattice3D<T,Descriptor>& lattice,
            Dot3D const& guoNode,
            std::vector<std::pair<int,int> > const& dryNodeFluidDirections,
            std::vector<plint> const& dryNodeIds, Dot3D const& absoluteOffset,
            const std::pair<int,int> &xDerivDirAndOrder, 
            const std::pair<int,int> &yDerivDirAndOrder, 
            const std::pair<int,int> &zDerivDirAndOrder,
            Array<T,3>& localForce, std::vector<AtomicBlock3D *> const& args );
private:
    bool useAllDirections;
public:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class GuoOffLatticeInfo3D : public ContainerBlockData {
    public:
        std::vector<Dot3D> const&                               getDryNodes() const
        { return dryNodes; }
        std::vector<Dot3D>&                                     getDryNodes()
        { return dryNodes; }
        std::vector<std::vector<std::pair<int,int> > > const&   getDryNodeFluidDirections() const
        { return dryNodeFluidDirections; }
        std::vector<std::vector<std::pair<int,int> > >&         getDryNodeFluidDirections()
        { return dryNodeFluidDirections; }
        std::vector<std::vector<plint> > const&                 getDryNodeIds() const
        { return dryNodeIds; }
        std::vector<std::vector<plint> >&                       getDryNodeIds()
        { return dryNodeIds; }
        std::vector<std::pair<int,int> > const&                 getXderivDirAndOrder() const
        { return xDerivDirAndOrder; }
        std::vector<std::pair<int,int> >&                       getXderivDirAndOrder()
        { return xDerivDirAndOrder; }
        std::vector<std::pair<int,int> > const&                 getYderivDirAndOrder() const
        { return yDerivDirAndOrder; }
        std::vector<std::pair<int,int> >&                       getYderivDirAndOrder()
        { return yDerivDirAndOrder; }
        std::vector<std::pair<int,int> > const&                 getZderivDirAndOrder() const
        { return zDerivDirAndOrder; }
        std::vector<std::pair<int,int> >&                       getZderivDirAndOrder()
        { return zDerivDirAndOrder; }
        Array<T,3> const&                                       getLocalForce() const
        { return localForce; }
        Array<T,3>&                                             getLocalForce()
        { return localForce; }
        virtual GuoOffLatticeInfo3D* clone() const {
            return new GuoOffLatticeInfo3D(*this);
        }
    private:
        std::vector<Dot3D> dryNodes;
        std::vector<std::vector<std::pair<int,int> > >   dryNodeFluidDirections;
        std::vector<std::vector<plint> >                 dryNodeIds;
        std::vector<std::pair<int,int> > xDerivDirAndOrder, yDerivDirAndOrder, zDerivDirAndOrder;
        Array<T,3>                                       localForce;
    };
};

}  // namespace plb

#endif  // GUO_OFF_LATTICE_MODEL_3D_H

