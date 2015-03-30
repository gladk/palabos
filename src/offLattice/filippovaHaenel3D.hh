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

#ifndef FILIPPOVA_HAENEL_3D_HH
#define FILIPPOVA_HAENEL_3D_HH

#include "offLattice/filippovaHaenel3D.h"
#include "offLattice/nextNeighbors3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "core/dynamics.h"
#include <algorithm>
#include <vector>
#include <cmath>

namespace plb {

template<typename T, template<typename U> class Descriptor>
FilippovaHaenelModel3D<T,Descriptor>::FilippovaHaenelModel3D (
        BoundaryShape3D<T,Array<T,3> >* shape_, int flowType_, bool useAllDirections_ )
    : OffLatticeModel3D<T,Array<T,3> >(shape_, flowType_),
      computeStat(true)
{ }

template<typename T, template<typename U> class Descriptor>
FilippovaHaenelModel3D<T,Descriptor>* FilippovaHaenelModel3D<T,Descriptor>::clone() const {
    return new FilippovaHaenelModel3D(*this);
}

template<typename T, template<typename U> class Descriptor>
plint FilippovaHaenelModel3D<T,Descriptor>::getNumNeighbors() const {
    return 1;
}

template<typename T, template<typename U> class Descriptor>
void FilippovaHaenelModel3D<T,Descriptor>::prepareCell (
        Dot3D const& cellLocation,
        AtomicContainerBlock3D& container )
{
    typedef Descriptor<T> D;
    Dot3D offset = container.getLocation();
    OffLatticeInfo3D* info = dynamic_cast<OffLatticeInfo3D*>(container.getData());
    PLB_ASSERT( info );
    std::vector<int> liquidNeighbors;
    std::vector<plint> ids;
    if (!this->isFluid(cellLocation+offset)) {
        for (int iPop=0; iPop<D::q; ++iPop) {
            Dot3D neighbor(cellLocation.x+D::c[iPop][0], cellLocation.y+D::c[iPop][1], cellLocation.z+D::c[iPop][2]);
            // If the non-fluid node has a fluid neighbor ...
            if (this->isFluid(neighbor+offset)) {
                // ... check how many fluid nodes it has ahead of it ...
                plint iTriangle=-1;
                global::timer("intersect").start();
                Array<T,3> locatedPoint;
                T distance;
                Array<T,3> wallNormal;
                Array<T,3> surfaceData;
                OffBoundary::Type bdType;
#ifdef PLB_DEBUG
                bool ok =
#endif
                    this->pointOnSurface (
                            cellLocation+offset, Dot3D(D::c[iPop][0],D::c[iPop][1],D::c[iPop][2]), locatedPoint, distance,
                            wallNormal, surfaceData, bdType, iTriangle );
                // In the following, the importance of directions is sorted wrt. how well they
                //   are aligned with the wall normal. It is better to take the continuous normal,
                //   because it is not sensitive to the choice of the triangle when we shoot at
                //   an edge.
                //wallNormal = this->computeContinuousNormal(locatedPoint, iTriangle);
                global::timer("intersect").stop();
                PLB_ASSERT( ok );
                // ... then add this node to the list.
                liquidNeighbors.push_back(iPop);
                ids.push_back(iTriangle);
            }
        }
        if (!liquidNeighbors.empty()) {
            info->getDryNodes().push_back(cellLocation);
            info->getDryNodeFluidDirections().push_back(liquidNeighbors);
            info->getDryNodeIds().push_back(ids);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ContainerBlockData*
    FilippovaHaenelModel3D<T,Descriptor>::generateOffLatticeInfo() const
{
    return new OffLatticeInfo3D;
}

template<typename T, template<typename U> class Descriptor>
Array<T,3> FilippovaHaenelModel3D<T,Descriptor>::getLocalForce (
                AtomicContainerBlock3D& container ) const
{
    OffLatticeInfo3D* info =
        dynamic_cast<OffLatticeInfo3D*>(container.getData());
    PLB_ASSERT( info );
    return info->getLocalForce();
}

template<typename T, template<typename U> class Descriptor>
void FilippovaHaenelModel3D<T,Descriptor>::boundaryCompletion (
        AtomicBlock3D& nonTypeLattice,
        AtomicContainerBlock3D& container,
        std::vector<AtomicBlock3D const*> const& args )
{
    BlockLattice3D<T,Descriptor>& lattice =
        dynamic_cast<BlockLattice3D<T,Descriptor>&> (nonTypeLattice);
    OffLatticeInfo3D* info =
        dynamic_cast<OffLatticeInfo3D*>(container.getData());
    PLB_ASSERT( info );
    std::vector<Dot3D> const&
        dryNodes = info->getDryNodes();
    std::vector<std::vector<int > > const&
        dryNodeFluidDirections = info->getDryNodeFluidDirections();
    std::vector<std::vector<plint> > const&
        dryNodeIds = info->getDryNodeIds();
    PLB_ASSERT( dryNodes.size() == dryNodeFluidDirections.size() );

    Dot3D absoluteOffset = lattice.getLocation();

    Array<T,3>& localForce = info->getLocalForce();
    localForce.resetToZero();
    for (pluint iDry=0; iDry<dryNodes.size(); ++iDry) {
        cellCompletion (
            lattice, dryNodes[iDry], dryNodeFluidDirections[iDry],
            dryNodeIds[iDry], absoluteOffset, localForce, args );
    }
}

template<typename T, template<typename U> class Descriptor>
void FilippovaHaenelModel3D<T,Descriptor>::cellCompletion (
        BlockLattice3D<T,Descriptor>& lattice, Dot3D const& guoNode,
        std::vector<int> const& dryNodeFluidDirections,
        std::vector<plint> const& dryNodeIds, Dot3D const& absoluteOffset,
        Array<T,3>& localForce, std::vector<AtomicBlock3D const*> const& args )
{
    typedef Descriptor<T> D;
    Cell<T,Descriptor>& s_cell =
        lattice.get( guoNode.x, guoNode.y, guoNode.z );
#ifdef PLB_DEBUG
    int noDynId =
#endif
        NoDynamics<T,Descriptor>().getId();
    PLB_ASSERT( s_cell.getDynamics().getId() == noDynId );
    for (plint iDirection=0; iDirection<(plint)dryNodeFluidDirections.size(); ++iDirection)
    {
        int iOpp = dryNodeFluidDirections[iDirection];
        int iPop = indexTemplates::opposite<Descriptor<T> >(iOpp);
        Dot3D fluidDirection(D::c[iOpp][0],D::c[iOpp][1],D::c[iOpp][2]);
        plint dryNodeId = dryNodeIds[iDirection];

        Array<T,3> wallNode, wall_vel;
        T wallDistance;
        OffBoundary::Type bdType;

        Cell<T,Descriptor> const& f_cell =
            lattice.get( guoNode.x+fluidDirection.x,
                         guoNode.y+fluidDirection.y,
                         guoNode.z+fluidDirection.z );

        Cell<T,Descriptor> const& ff_cell =
            lattice.get( guoNode.x+2*fluidDirection.x,
                         guoNode.y+2*fluidDirection.y,
                         guoNode.z+2*fluidDirection.z );

        Cell<T,Descriptor> collidedCell(f_cell);
        BlockStatistics statsCopy(lattice.getInternalStatistics());
        collidedCell.collide(statsCopy);

        T f_rhoBar, ff_rhoBar;
        Array<T,3> f_j, ff_j;
        Array<T,3> wallNormal;
        f_cell.getDynamics().computeRhoBarJ(f_cell, f_rhoBar, f_j);
        ff_cell.getDynamics().computeRhoBarJ(ff_cell, ff_rhoBar, ff_j);
        T f_rho = D::fullRho(f_rhoBar);
        T f_jSqr = normSqr(f_j);

#ifdef PLB_DEBUG
        bool ok =
#endif
        this->pointOnSurface( guoNode+absoluteOffset, fluidDirection,
                              wallNode, wallDistance, wallNormal,
                              wall_vel, bdType, dryNodeId );
        PLB_ASSERT( ok );

        Array<T,3> w_j = wall_vel*f_rho;
        T d = std::sqrt(D::cNormSqr[iOpp]);
        PLB_ASSERT( wallDistance <= d );
        T delta = 1.0-wallDistance / d;

        T kappa = 0.;
        Array<T,3> wf_j; wf_j.resetToZero();
        T omega = f_cell.getDynamics().getOmega();


        if (delta<0.5) {
            //wf_j = f_j;
            //kappa = (omega*(2.0*delta-1.0))/(1.0-omega);

            wf_j = ff_j;
            kappa = (omega*(2*delta-1.0))/(1.0-2.0*omega);
        }
        else {
            //wf_j = f_j * ((delta-1.0)/delta) + w_j / delta;
            //kappa = omega*(2.0*delta-1.0)

            wf_j = (1.0-3.0/(2.0*delta))*f_j+3.0/(2.0*delta)*w_j;
            kappa = (2.0*omega*(2.0*delta-1.0))/(2.0+omega);
        }

        T c_i_wf_j_f_j = D::c[iPop][0]*(wf_j[0]-f_j[0]) +
                         D::c[iPop][1]*(wf_j[1]-f_j[1]) +
                         D::c[iPop][2]*(wf_j[2]-f_j[2]) ;

        T c_i_w_j     = D::c[iPop][0]*w_j[0] + D::c[iPop][1]*w_j[1] + D::c[iPop][2]*w_j[2];


        T f_ieq = f_cell.getDynamics().computeEquilibrium(iPop, f_rhoBar, f_j, f_jSqr) +
                      D::t[iPop]*D::invCs2*c_i_wf_j_f_j;
        s_cell[iOpp] = (1.0-kappa)*collidedCell[iPop]+kappa*f_ieq+2.0*D::t[iPop]*D::invCs2*c_i_w_j;

        localForce[0] += D::c[iPop][0]*(s_cell[iPop]+s_cell[iOpp]);
        localForce[1] += D::c[iPop][1]*(s_cell[iPop]+s_cell[iOpp]);
        localForce[2] += D::c[iPop][2]*(s_cell[iPop]+s_cell[iOpp]);
    }
}

}  // namespace plb

#endif  // FILIPPOVA_HAENEL_3D_HH

