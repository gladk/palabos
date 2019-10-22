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
 * Groups all the include files for the 3D off-lattice directory.
 */

#include "offLattice/marchingCube.h"
#include "offLattice/triangleSelector.h"
#include "offLattice/triangleSet.h"
#include "offLattice/connectedTriangleSet.h"
#include "offLattice/nextNeighbors3D.h"
#include "offLattice/triangleToDef.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "offLattice/voxelizer.h"
#include "offLattice/makeSparse3D.h"
#include "offLattice/triangleHash.h"
#include "offLattice/offLatticeBoundaryProcessor3D.h"
#include "offLattice/offLatticeBoundaryProfiles3D.h"
#include "offLattice/offLatticeBoundaryCondition3D.h"
#include "offLattice/boundaryShapes3D.h"
#include "offLattice/triangleBoundary3D.h"
#include "offLattice/offLatticeModel3D.h"
#include "offLattice/guoOffLatticeModel3D.h"
#include "offLattice/bouzidiOffLatticeModel3D.h"
#include "offLattice/guoAdvDiffOffLatticeModel3D.h"
#include "offLattice/triangleSetGenerator.h"
#include "offLattice/immersedWalls3D.h"
#include "offLattice/immersedAdvectionDiffusionWalls3D.h"
#include "offLattice/filippovaHaenelOffLatticeModel3D.h"

#ifndef PLB_BGP
#ifdef PLB_USE_EIGEN
#include "offLattice/generalizedOffLatticeModel3D.h"
#endif
#endif

