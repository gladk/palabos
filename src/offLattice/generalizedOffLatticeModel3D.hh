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

#ifndef GENERALIZED_OFF_LATTICE_MODEL_3D_HH
#define GENERALIZED_OFF_LATTICE_MODEL_3D_HH

#include "offLattice/generalizedOffLatticeModel3D.h"
#include "offLattice/nextNeighbors3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include <algorithm>
#include <cmath>

namespace plb {

template<typename T, template<typename U> class Descriptor>
ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>::ExtrapolatedGeneralizedOffLatticeModel3D (
        BoundaryShape3D<T,Array<T,3> >* shape_, int flowType_ )
    : OffLatticeModel3D<T,Array<T,3> >(shape_, flowType_)
{ }

template<typename T, template<typename U> class Descriptor>
ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>::ExtrapolatedGeneralizedOffLatticeModel3D (
        ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor> const& rhs )
    : OffLatticeModel3D<T,Array<T,3> >(rhs)
{ }

template<typename T, template<typename U> class Descriptor>
ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>&
    ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>::operator=(ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor> const& rhs)
{
    OffLatticeModel3D<T,Array<T,3> >::operator=(rhs);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>* ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>::clone() const {
    return new ExtrapolatedGeneralizedOffLatticeModel3D(*this);
}

template<typename T, template<typename U> class Descriptor>
plint ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>::getNumNeighbors() const {
    return 2;
}

template<typename T, template<typename U> class Descriptor>
void ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>::prepareCell (
        Dot3D const& cellLocation,
        AtomicContainerBlock3D& container )
{
    Dot3D offset = container.getLocation();
    ExtrapolatedGeneralizedOffLatticeInfo3D* info =
        dynamic_cast<ExtrapolatedGeneralizedOffLatticeInfo3D*>(container.getData());
    PLB_ASSERT( info );
    if (!this->isFluid(cellLocation+offset)) {
        std::vector<std::pair<int,int> > liquidNeighbors;
        std::vector<int> liquidNeighborsNoSolid;
        std::vector<plint> ids;
        for (int iNeighbor=0; iNeighbor<NextNeighbor<T>::numNeighbors; ++iNeighbor) {
            int const* c = NextNeighbor<T>::c[iNeighbor];
            Dot3D neighbor(cellLocation.x+c[0], cellLocation.y+c[1], cellLocation.z+c[2]);
            // If the non-fluid node has a fluid neighbor ...
            if (this->isFluid(neighbor+offset)) {
                // ... check how many fluid nodes it has ahead of it ...
                int depth = 1;
                for (int iDepth=2; iDepth<=getNumNeighbors(); ++iDepth) {
                    Dot3D nextNeighbor(cellLocation.x+iDepth*c[0],
                                       cellLocation.y+iDepth*c[1],
                                       cellLocation.z+iDepth*c[2]);
                    if (this->isFluid(nextNeighbor+offset)) {
                        depth = iDepth;
                    }
                    else {
                        break;
                    }
                }
                // ... then add this node to the list.
                liquidNeighbors.push_back(std::make_pair(iNeighbor,depth));
                Array<T,3> locatedPoint;
                T distance;
                Array<T,3> wallNormal;
                Array<T,3> surfaceData;
                plint iTriangle=-1;
                OffBoundary::Type bdType;
#ifdef PLB_DEBUG
                bool ok =
#endif
                    this->pointOnSurface (
                            cellLocation+offset, Dot3D(c[0],c[1],c[2]), locatedPoint, distance,
                            wallNormal, surfaceData, bdType, iTriangle );
                PLB_ASSERT( ok );
                ids.push_back(iTriangle);
                liquidNeighborsNoSolid.push_back(iNeighbor);
            }
            else {
                bool fluidNeighbor = false;
                for (int jNeighbor=1; jNeighbor<Descriptor<T>::q; ++jNeighbor) {
                    Dot3D next_neighbor(neighbor.x+Descriptor<T>::c[jNeighbor][0], 
                                        neighbor.y+Descriptor<T>::c[jNeighbor][1], 
                                        neighbor.z+Descriptor<T>::c[jNeighbor][2]);
                    if ((this->isFluid(next_neighbor+offset))) {
                        fluidNeighbor = true;
                        break;
                    }
                }
                if (fluidNeighbor) {
                    liquidNeighborsNoSolid.push_back(iNeighbor);
                }
            }
        }
        if (!liquidNeighbors.empty()) {
            // selecting only deepest depth available (for a higher order interpolation order)
//             std::vector<std::vector<std::pair<int,int> > > selectionLiquid(getNumNeighbors()+1);
//             std::vector<std::vector<plint > > selectionIds(getNumNeighbors()+1);
//             for (pluint iA = 0; iA < liquidNeighbors.size(); ++iA) {
//                 selectionLiquid[liquidNeighbors[iA].second].push_back(liquidNeighbors[iA]);
//                 selectionIds[liquidNeighbors[iA].second].push_back(ids[iA]);
//             }
            
            info->getDryNodes().push_back(cellLocation);
            info->getDryNodeFluidWithFluidDirections().push_back(liquidNeighborsNoSolid);
            info->getDryNodeFluidDirections().push_back(liquidNeighbors);
            info->getDryNodeIds().push_back(ids);
            
            // add only biggest depth to the list of directions and triangles
//             for (plint iA = getNumNeighbors(); iA >= 0; --iA) {
//                 if (!selectionLiquid[iA].empty()) {
//                     info->getDryNodeFluidDirections().push_back(selectionLiquid[iA]);
//                     info->getDryNodeIds().push_back(selectionIds[iA]);
//                     break;
//                 }
//             }
            
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ContainerBlockData*
    ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>::generateOffLatticeInfo() const
{
    return new ExtrapolatedGeneralizedOffLatticeInfo3D;
}

template<typename T, template<typename U> class Descriptor>
Array<T,3> ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>::getLocalForce (
                AtomicContainerBlock3D& container ) const
{
    ExtrapolatedGeneralizedOffLatticeInfo3D* info =
        dynamic_cast<ExtrapolatedGeneralizedOffLatticeInfo3D*>(container.getData());
    PLB_ASSERT( info );
    return info->getLocalForce();
}

template<typename T, template<typename U> class Descriptor>
void ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>::boundaryCompletion (
        AtomicBlock3D& nonTypeLattice,
        AtomicContainerBlock3D& container,
        std::vector<AtomicBlock3D const*> const& args )
{
    BlockLattice3D<T,Descriptor>& lattice =
        dynamic_cast<BlockLattice3D<T,Descriptor>&> (nonTypeLattice);
    ExtrapolatedGeneralizedOffLatticeInfo3D* info =
        dynamic_cast<ExtrapolatedGeneralizedOffLatticeInfo3D*>(container.getData());
    PLB_ASSERT( info );
    std::vector<Dot3D> const&
        dryNodes = info->getDryNodes();
    std::vector<std::vector<std::pair<int,int> > > const&
        dryNodeFluidDirections = info->getDryNodeFluidDirections();
    std::vector<std::vector<int> > const&
        dryNodeFluidNoSolidDirections = info->getDryNodeFluidWithFluidDirections();
    std::vector<std::vector<plint> > const&
        dryNodeIds = info->getDryNodeIds();
    PLB_ASSERT( dryNodes.size() == dryNodeFluidDirections.size() );

    Dot3D absoluteOffset = lattice.getLocation();

    Array<T,3>& localForce = info->getLocalForce();
    localForce.resetToZero();
    for (pluint iDry=0; iDry<dryNodes.size(); ++iDry) {
        cellCompletion (
            lattice, dryNodes[iDry], dryNodeFluidDirections[iDry], dryNodeFluidNoSolidDirections[iDry],
            dryNodeIds[iDry], absoluteOffset, localForce );
    }
}

template<typename T, template<typename U> class Descriptor>
void ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>::cellCompletion (
        BlockLattice3D<T,Descriptor>& lattice,
        Dot3D const& genNode,
        std::vector<std::pair<int,int> > const& dryNodeFluidDirections,
        std::vector<int> const& dryNodeFluidNoSolidDirections,
        std::vector<plint> const& dryNodeIds, Dot3D const& absoluteOffset,
        Array<T,3>& localForce )
{
    typedef Descriptor<T> D;
    using namespace indexTemplates;
    
    Cell<T,Descriptor>& cell =
        lattice.get(genNode.x, genNode.y, genNode.z);
    plint numDirections = (plint)dryNodeFluidDirections.size();
    std::vector<T> weights(numDirections);
    std::vector<Array<T,Descriptor<T>::d> > u_vect(numDirections);
    T sumWeights = T();
    Array<T,3> wallNormal;
    for (plint iDirection=0; iDirection<numDirections; ++iDirection) {
        int iNeighbor = dryNodeFluidDirections[iDirection].first;
        int const* c = NextNeighbor<T>::c[iNeighbor];
        Dot3D fluidDirection(c[0],c[1],c[2]);
        plint dryNodeId = dryNodeIds[iDirection];
        int depth = dryNodeFluidDirections[iDirection].second;

        Array<T,3> wallNode, wall_vel;
        T wallDistance;
        OffBoundary::Type bdType;
#ifdef PLB_DEBUG
        bool ok =
#endif
        this->pointOnSurface( genNode+absoluteOffset, fluidDirection,
                        wallNode, wallDistance, wallNormal,
                        wall_vel, bdType, dryNodeId );
        PLB_ASSERT( ok );
        for (int iD=0; iD<Descriptor<T>::d; ++iD) {
            // Use the formula uLB = uP - 1/2 g. If there is no external force,
            //   the force term automatically evaluates to zero.
            wall_vel[iD] -= (T)0.5*getExternalForceComponent(cell,iD);
        }
        T invDistanceToNeighbor = NextNeighbor<T>::invD[iNeighbor];
        PLB_ASSERT( wallDistance <= NextNeighbor<T>::d[iNeighbor] );
        T delta = (T)1. - wallDistance * invDistanceToNeighbor;
        Array<T,3> normalFluidDirection((T)fluidDirection.x, (T)fluidDirection.y, (T)fluidDirection.z);
        normalFluidDirection *= invDistanceToNeighbor;
        weights[iDirection] = std::fabs ( dot(normalFluidDirection, wallNormal) );
        sumWeights += weights[iDirection];

        computeVelocity(
                lattice, genNode, fluidDirection, depth,
                wallNode, delta, wall_vel, bdType, wallNormal,
                u_vect[iDirection] );
    }

    Array<T,D::d> u; u.resetToZero();
    for (plint iDirection=0; iDirection<numDirections; ++iDirection) {
        u += u_vect[iDirection] * weights[iDirection];
    }
    u /= sumWeights;
    
    Array<T,D::d> deltaJ;
    deltaJ.resetToZero();
    for (plint iDirection=0; iDirection<numDirections; ++iDirection) {
        int iNeighbor = dryNodeFluidDirections[iDirection].first;
        int iPop = nextNeighborPop<T,Descriptor>(iNeighbor);
        if (iPop>=0) {
            plint oppPop = indexTemplates::opposite<D>(iPop);
            deltaJ[0] += D::c[oppPop][0]*cell[oppPop];
            deltaJ[1] += D::c[oppPop][1]*cell[oppPop];
            deltaJ[2] += D::c[oppPop][2]*cell[oppPop];
        }
    }
    
    std::vector<plint> knownIndices;
    knownIndices.push_back(0);
    for (pluint iDirection=0; iDirection<dryNodeFluidNoSolidDirections.size(); ++iDirection) {
        int iNeighbor = dryNodeFluidNoSolidDirections[iDirection];
        plint iPop = nextNeighborPop<T,Descriptor>(iNeighbor);
        if (iPop>=0) {
            plint index = opposite<Descriptor<T> >(iPop);
            knownIndices.push_back(index);
        }
    }
    PLB_ASSERT(knownIndices.size() >= 6);
    
    std::vector<plint> missingIndices = remainingIndexes<Descriptor<T> >(knownIndices);
    
    Dynamics<T,Descriptor> const& dynamics = cell.getDynamics();
    DirichletVelocityBoundarySolver<T,Descriptor> bc(missingIndices, knownIndices, u);
    bc.apply(cell,dynamics,!this->getPartialReplace());
    
    Cell<T,Descriptor> cellCopy(cell);
    BlockStatistics statsCopy(lattice.getInternalStatistics());
    cellCopy.collide(statsCopy);
    
    for (plint iDirection=0; iDirection<(plint)dryNodeFluidDirections.size(); ++iDirection) {
        int iNeighbor = dryNodeFluidDirections[iDirection].first;
        plint iPop = nextNeighborPop<T,Descriptor>(iNeighbor);
        if (iPop>=0) {
            deltaJ[0] -= D::c[iPop][0]*cellCopy[iPop];
            deltaJ[1] -= D::c[iPop][1]*cellCopy[iPop];
            deltaJ[2] -= D::c[iPop][2]*cellCopy[iPop];
        }
    }
    localForce += deltaJ;
}

template<typename T, template<typename U> class Descriptor>
void ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor>::computeVelocity (
          BlockLattice3D<T,Descriptor> const& lattice, Dot3D const& genNode,
          Dot3D const& fluidDirection, int depth, Array<T,3> const& wallNode, T delta,
          Array<T,3> const& wall_u, OffBoundary::Type bdType, Array<T,3> const& wallNormal,
          Array<T,Descriptor<T>::d>& u) const
{
    Array<T,Descriptor<T>::d> u1, u2;
    Cell<T,Descriptor> const& cell1 =
        lattice.get( genNode.x+fluidDirection.x,
                     genNode.y+fluidDirection.y,
                     genNode.z+fluidDirection.z );
    Cell<T,Descriptor> const& cell2 =
        lattice.get( genNode.x+2*fluidDirection.x,
                     genNode.y+2*fluidDirection.y,
                     genNode.z+2*fluidDirection.z );
    if (!this->velIsJ()) {
        cell1.getDynamics().computeVelocity(cell1, u1);
        cell2.getDynamics().computeVelocity(cell2, u2);
    } else {
        T rhoBar;
        cell1.getDynamics().computeRhoBarJ(cell1, rhoBar, u1);
        cell2.getDynamics().computeRhoBarJ(cell2, rhoBar, u2);
    }

    //bool neumann = bdType==OffBoundary::neumann;
    if (depth < 2) {
        if (delta < (T)0.25) {
            u = wall_u;
        }
        else {
            u = (T)1./delta * (wall_u+(delta-(T)1.)*u1);
        }
    }
    else {  // depth >= 2
        if (delta < (T)0.75) {
            u = wall_u + (delta-(T)1.)*u1 +
                ((T)1.-delta)/((T)1.+delta)*((T)2.*wall_u+(delta-(T)1.)*u2);
        }
        else {
            u = (T)1./delta * (wall_u+(delta-(T)1.)*u1);
        }
    }
    if ( bdType==OffBoundary::neumann )
    {
        u = u1;
    }
    if ( bdType==OffBoundary::densityNeumann ) {
        u = dot(u1,wallNormal)*wallNormal;
    }
}


// ========================================================================= //
// ====Generalized BC with the velocity ocmputed with an interpolation====== //
// ========================================================================= //


template<typename T, template<typename U> class Descriptor>
InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>::InterpolatedGeneralizedOffLatticeModel3D (
        BoundaryShape3D<T,Array<T,3> >* shape_, int flowType_ )
    : OffLatticeModel3D<T,Array<T,3> >(shape_, flowType_)
{ }

template<typename T, template<typename U> class Descriptor>
InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>::InterpolatedGeneralizedOffLatticeModel3D (
        InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor> const& rhs )
    : OffLatticeModel3D<T,Array<T,3> >(rhs)
{ }

template<typename T, template<typename U> class Descriptor>
InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>&
    InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>::operator=(InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor> const& rhs)
{
    OffLatticeModel3D<T,Array<T,3> >::operator=(rhs);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>* InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>::clone() const {
    return new InterpolatedGeneralizedOffLatticeModel3D(*this);
}

template<typename T, template<typename U> class Descriptor>
plint InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>::getNumNeighbors() const {
    return 2;
}

template<typename T, template<typename U> class Descriptor>
void InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>::prepareCell (
        Dot3D const& cellLocation,
        AtomicContainerBlock3D& container )
{
    Dot3D offset = container.getLocation();
    InterpolatedGeneralizedOffLatticeInfo3D* info =
        dynamic_cast<InterpolatedGeneralizedOffLatticeInfo3D*>(container.getData());
    PLB_ASSERT( info );
    if (this->isFluid(cellLocation+offset)) {
        std::vector<std::pair<int,int> > solidNeighbors;
        std::vector<int> wetNodeFluidDirections;
        std::vector<plint> ids;
        for (int iNeighbor=0; iNeighbor<NextNeighbor<T>::numNeighbors; ++iNeighbor) {
            int const* c = NextNeighbor<T>::c[iNeighbor];
            Dot3D neighbor(cellLocation.x+c[0], cellLocation.y+c[1], cellLocation.z+c[2]);
            // TODO use only depth two if possible (add a check for only depth 
            // two selection or depth one if depth two not possible.
            
            // If the fluid node has a non-fluid neighbor ...
            if (!this->isFluid(neighbor+offset)) {
                // ... check how many fluid nodes without solid neighbors 
                // in the direction opposite to the solid.
                int depth = 0;
                for (int iDepth=1; iDepth<=getNumNeighbors(); ++iDepth) {
                    Dot3D nextNeighbor(cellLocation.x-iDepth*c[0],
                                       cellLocation.y-iDepth*c[1],
                                       cellLocation.z-iDepth*c[2]);
                    
                    bool solidNeighbor = false;
                    for (int iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                        Dot3D nextNeighbor_neighbor(nextNeighbor.x+Descriptor<T>::c[iPop][0],
                                                    nextNeighbor.y+Descriptor<T>::c[iPop][1],
                                                    nextNeighbor.z+Descriptor<T>::c[iPop][2]);

                        if (!(this->isFluid(nextNeighbor_neighbor+offset))) {
                            solidNeighbor = true;
                            break;
                        }
                    }
                    
                    if (this->isFluid(nextNeighbor+offset) && !solidNeighbor) {
                        depth = iDepth;
                    }
                    else {
                        break;
                    }
                }
                // ... then add this node to the list.
                solidNeighbors.push_back(std::make_pair(iNeighbor,depth));
                Array<T,3> locatedPoint;
                T distance;
                Array<T,3> wallNormal;
                Array<T,3> surfaceData;
                plint iTriangle=-1;
                OffBoundary::Type bdType;
#ifdef PLB_DEBUG
                bool ok =
#endif
                    this->pointOnSurface (
                            cellLocation+offset, Dot3D(c[0],c[1],c[2]), locatedPoint, distance,
                            wallNormal, surfaceData, bdType, iTriangle );
                PLB_ASSERT( ok );
                ids.push_back(iTriangle);
            }
        }
        
        if (!solidNeighbors.empty()) {
            // selecting only deepest depth available (for a higher order interpolation order)
            std::vector<std::vector<std::pair<int,int> > > selectionSolid(getNumNeighbors()+1);
            std::vector<std::vector<plint > > selectionIds(getNumNeighbors()+1);
            for (pluint iA = 0; iA < solidNeighbors.size(); ++iA) {
                selectionSolid[solidNeighbors[iA].second].push_back(solidNeighbors[iA]);
                selectionIds[solidNeighbors[iA].second].push_back(ids[iA]);
            }
            // adding fluid boundary wet node's fluid directions (know populations directions).
            for (int iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                int const* c = Descriptor<T>::c[iPop];
                Dot3D potFluidNeighbor(cellLocation.x+c[0], cellLocation.y+c[1], cellLocation.z+c[2]);
                
                if (this->isFluid(potFluidNeighbor+offset)) {
                    wetNodeFluidDirections.push_back(iPop);
                }
            }
            
            info->getWetNodes().push_back(cellLocation);
            // add only biggest depth to the list of directions and triangles
            for (plint iA = getNumNeighbors(); iA >= 0; --iA) {
                if (!selectionSolid[iA].empty()) {
                    info->getWetNodeSolidDirections().push_back(selectionSolid[iA]);
                    info->getWetNodeIds().push_back(selectionIds[iA]);
                    break;
                }
            }
            info->getWetNodeFluidDirections().push_back(wetNodeFluidDirections);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ContainerBlockData*
    InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>::generateOffLatticeInfo() const
{
    return new InterpolatedGeneralizedOffLatticeInfo3D;
}

template<typename T, template<typename U> class Descriptor>
Array<T,3> InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>::getLocalForce (
                AtomicContainerBlock3D& container ) const
{
    InterpolatedGeneralizedOffLatticeInfo3D* info =
        dynamic_cast<InterpolatedGeneralizedOffLatticeInfo3D*>(container.getData());
    PLB_ASSERT( info );
    return info->getLocalForce();
}

template<typename T, template<typename U> class Descriptor>
void InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>::boundaryCompletion (
        AtomicBlock3D& nonTypeLattice,
        AtomicContainerBlock3D& container,
        std::vector<AtomicBlock3D const*> const& args )
{
    BlockLattice3D<T,Descriptor>& lattice =
        dynamic_cast<BlockLattice3D<T,Descriptor>&> (nonTypeLattice);
    InterpolatedGeneralizedOffLatticeInfo3D* info =
        dynamic_cast<InterpolatedGeneralizedOffLatticeInfo3D*>(container.getData());
    PLB_ASSERT( info );
    std::vector<Dot3D> const&
        wetNodes = info->getWetNodes();
    std::vector<std::vector<std::pair<int,int> > > const&
        wetNodeSolidDirections = info->getWetNodeSolidDirections();
    std::vector<std::vector<int> > const&
        wetNodeFluidDirections = info->getWetNodeFluidDirections();
    std::vector<std::vector<plint> > const&
        wetNodeIds = info->getWetNodeIds();
    PLB_ASSERT( wetNodes.size() == wetNodeSolidDirections.size() );

    Dot3D absoluteOffset = lattice.getLocation();

    Array<T,3>& localForce = info->getLocalForce();
    localForce.resetToZero();
    for (pluint iWet=0; iWet<wetNodes.size(); ++iWet) {
        cellCompletion (
            lattice, wetNodes[iWet], wetNodeSolidDirections[iWet], wetNodeFluidDirections[iWet],
            wetNodeIds[iWet], absoluteOffset, localForce );
    }
}

template<typename T, template<typename U> class Descriptor>
void InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>::cellCompletion (
        BlockLattice3D<T,Descriptor>& lattice,
        Dot3D const& genNode,
        std::vector<std::pair<int,int> > const& wetNodeSolidDirections,
        std::vector<int> const& wetNodeFluidDirections,
        std::vector<plint> const& wetNodeIds, Dot3D const& absoluteOffset,
        Array<T,3>& localForce )
{
    typedef Descriptor<T> D;
    using namespace indexTemplates;
    
    Cell<T,Descriptor>& cell =
        lattice.get(genNode.x, genNode.y, genNode.z);
    plint numDirections = (plint)wetNodeSolidDirections.size();
    std::vector<T> weights(numDirections);
    std::vector<Array<T,Descriptor<T>::d> > u_vect(numDirections);
    T sumWeights = T();
    Array<T,3> wallNormal;
    for (plint iDirection=0; iDirection<numDirections; ++iDirection) {
        int iNeighbor = wetNodeSolidDirections[iDirection].first;
        int depth = wetNodeSolidDirections[iDirection].second;
//         pcout << "depth = " << depth << std::endl;
        int const* c = NextNeighbor<T>::c[iNeighbor];
        
        Dot3D solidDirection(c[0],c[1],c[2]);
        plint wetNodeId = wetNodeIds[iDirection];
        
        Array<T,3> wallNode, wall_vel;
        T wallDistance;
        OffBoundary::Type bdType;
#ifdef PLB_DEBUG
        bool ok =
#endif
        this->pointOnSurface( genNode+absoluteOffset, solidDirection,
                        wallNode, wallDistance, wallNormal,
                        wall_vel, bdType, wetNodeId );
        PLB_ASSERT( ok );
        for (int iD=0; iD<Descriptor<T>::d; ++iD) {
            // Use the formula uLB = uP - 1/2 g. If there is no external force,
            //   the force term automatically evaluates to zero.
            wall_vel[iD] -= (T)0.5*getExternalForceComponent(cell,iD);
        }
        T invDistanceToNeighbor = NextNeighbor<T>::invD[iNeighbor];
        PLB_ASSERT( wallDistance <= NextNeighbor<T>::d[iNeighbor] );
        
        Array<T,3> normalFluidDirection((T)solidDirection.x, (T)solidDirection.y, (T)solidDirection.z);
        normalFluidDirection *= invDistanceToNeighbor;
        weights[iDirection] = std::fabs ( dot(normalFluidDirection, wallNormal) );
        sumWeights += weights[iDirection];

        computeVelocity (
                lattice, genNode, solidDirection, depth,
                wallNode, wallDistance, NextNeighbor<T>::d[iNeighbor], wall_vel, bdType, wallNormal,
                u_vect[iDirection] );
    }

    Array<T,D::d> u; u.resetToZero();
    for (plint iDirection=0; iDirection<numDirections; ++iDirection) {
        u += u_vect[iDirection] * weights[iDirection];
    }
    u /= sumWeights;
    
//     pcout << "vel = " << u[0] << ", " << u[1] << ", " << u[2] << std::endl;

    Array<T,D::d> deltaJ;
    deltaJ.resetToZero();
    for (plint iDirection=0; iDirection<numDirections; ++iDirection) {
        int iNeighbor = wetNodeSolidDirections[iDirection].first;
        int iPop = nextNeighborPop<T,Descriptor>(iNeighbor);
        if (iPop>=0) {
            deltaJ[0] += D::c[iPop][0]*cell[iPop];
            deltaJ[1] += D::c[iPop][1]*cell[iPop];
            deltaJ[2] += D::c[iPop][2]*cell[iPop];
        }
    }
    
    std::vector<plint> knownIndices;
    knownIndices.push_back(0);
//     pcout << "indexes = " << std::endl;
    for (pluint iDirection=0; iDirection<wetNodeFluidDirections.size(); ++iDirection) {
        int iPop = wetNodeFluidDirections[iDirection];
        plint index = opposite<Descriptor<T> >(iPop);
        knownIndices.push_back(index);
//             pcout << index << std::endl;
    }
    PLB_ASSERT(knownIndices.size() >= 6);
    
    std::vector<plint> missingIndices = remainingIndexes<Descriptor<T> >(knownIndices);
    
    Dynamics<T,Descriptor> const& dynamics = cell.getDynamics();
    DirichletVelocityBoundarySolver<T,Descriptor> bc(missingIndices, knownIndices, u);
    bc.apply(cell,dynamics,!this->getPartialReplace());
    
    Cell<T,Descriptor> cellCopy(cell);
    BlockStatistics statsCopy(lattice.getInternalStatistics());
    cellCopy.collide(statsCopy);
    
    for (plint iDirection=0; iDirection<numDirections; ++iDirection) {
        int iNeighbor = wetNodeSolidDirections[iDirection].first;
        plint iPop = nextNeighborPop<T,Descriptor>(iNeighbor);
        if (iPop>=0) {
            plint oppPop = indexTemplates::opposite<D>(iPop);
            deltaJ[0] -= D::c[oppPop][0]*cellCopy[oppPop];
            deltaJ[1] -= D::c[oppPop][1]*cellCopy[oppPop];
            deltaJ[2] -= D::c[oppPop][2]*cellCopy[oppPop];
        }
    }
    localForce += deltaJ;
}

template<typename T, template<typename U> class Descriptor>
void InterpolatedGeneralizedOffLatticeModel3D<T,Descriptor>::computeVelocity (
          BlockLattice3D<T,Descriptor> const& lattice, Dot3D const& genNode,
          Dot3D const& solidDirection, int depth, Array<T,3> const& wallNode, T wallDistance, T cellDistance,
          Array<T,3> const& wall_u, OffBoundary::Type bdType, Array<T,3> const& wallNormal,
          Array<T,Descriptor<T>::d>& u) const
{
    //bool neumann = bdType==OffBoundary::neumann;
    Array<T,Descriptor<T>::d> u1, u2;
    if (depth == 0) {
        u = wall_u;
    }
    else if (depth == 1) {
        Cell<T,Descriptor> const& cell1 =
        lattice.get( genNode.x-solidDirection.x,
                     genNode.y-solidDirection.y,
                     genNode.z-solidDirection.z );
        
        if (!this->velIsJ()) {
            cell1.getDynamics().computeVelocity(cell1, u1);
        } else {
            T rhoBar;
            cell1.getDynamics().computeRhoBarJ(cell1, rhoBar, u1);
        }
        
        u = (wall_u * cellDistance + wallDistance * u1) / (wallDistance+cellDistance);
    }
    else if (depth >= 2) {  // depth >= 2
        Cell<T,Descriptor> const& cell1 =
        lattice.get( genNode.x-solidDirection.x,
                     genNode.y-solidDirection.y,
                     genNode.z-solidDirection.z );
        Cell<T,Descriptor> const& cell2 =
        lattice.get( genNode.x-2*solidDirection.x,
                     genNode.y-2*solidDirection.y,
                     genNode.z-2*solidDirection.z );
        if (!this->velIsJ()) {
            cell1.getDynamics().computeVelocity(cell1, u1);
            cell2.getDynamics().computeVelocity(cell2, u2);
        } else {
            T rhoBar;
            cell1.getDynamics().computeRhoBarJ(cell1, rhoBar, u1);
            cell2.getDynamics().computeRhoBarJ(cell2, rhoBar, u2);
        }
        
        T invDenom = (T)1 / (wallDistance*wallDistance+(T)2*cellDistance*cellDistance+(T)3*wallDistance*cellDistance);
        
        u = wallDistance*(-wallDistance-cellDistance)*u2*invDenom 
            + (T)2*wallDistance*(wallDistance+(T)2*cellDistance)*u1*invDenom
            + (T)2*cellDistance*cellDistance*wall_u*invDenom;
    }
    if ( bdType==OffBoundary::neumann )
    {
        Cell<T,Descriptor> const& cell1 =
        lattice.get( genNode.x-solidDirection.x,
                     genNode.y-solidDirection.y,
                     genNode.z-solidDirection.z );
        
        cell1.getDynamics().computeVelocity(cell1, u1);
        
        u = u1;
    }
    if ( bdType==OffBoundary::densityNeumann ) {
        Cell<T,Descriptor> const& cell1 =
        lattice.get( genNode.x-solidDirection.x,
                     genNode.y-solidDirection.y,
                     genNode.z-solidDirection.z );
        
        cell1.getDynamics().computeVelocity(cell1, u1);
        
        u = dot(u1,wallNormal)*wallNormal;
    }
    if ( bdType==OffBoundary::freeSlip ) {
        Cell<T,Descriptor> const& cell1 =
        lattice.get( genNode.x-solidDirection.x,
                     genNode.y-solidDirection.y,
                     genNode.z-solidDirection.z );
        
        cell1.getDynamics().computeVelocity(cell1, u1);
        
        u = u1-dot(u1,wallNormal)*wallNormal;
    }
}

}  // namespace plb

#endif  // GUO_OFF_LATTICE_MODEL_3D_HH
