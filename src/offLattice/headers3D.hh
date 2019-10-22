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
 * Groups all the template files for the 3D off-lattice directory.
 */

#include "offLattice/marchingCube.hh"
#include "offLattice/triangleSet.hh"
#include "offLattice/connectedTriangleSet.hh"
#include "offLattice/nextNeighbors3D.hh"
#include "offLattice/triangleToDef.hh"
#include "offLattice/triangularSurfaceMesh.hh"
#include "offLattice/voxelizer.hh"
#include "offLattice/makeSparse3D.hh"
#include "offLattice/triangleHash.hh"
#include "offLattice/offLatticeBoundaryProcessor3D.hh"
#include "offLattice/offLatticeBoundaryProfiles3D.hh"
#include "offLattice/offLatticeBoundaryCondition3D.hh"
#include "offLattice/boundaryShapes3D.hh"
#include "offLattice/triangleBoundary3D.hh"
#include "offLattice/offLatticeModel3D.hh"
#include "offLattice/guoOffLatticeModel3D.hh"
#include "offLattice/bouzidiOffLatticeModel3D.hh"
#include "offLattice/guoAdvDiffOffLatticeModel3D.hh"
#include "offLattice/triangleSetGenerator.hh"
#include "offLattice/immersedWalls3D.hh"
#include "offLattice/immersedAdvectionDiffusionWalls3D.hh"
#include "offLattice/filippovaHaenelOffLatticeModel3D.hh"

#ifndef PLB_BGP
#ifdef PLB_USE_EIGEN
#include "offLattice/generalizedOffLatticeModel3D.hh"
#endif
#endif

