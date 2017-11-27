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
 * Groups all the 3D include files in the directory multiPhysics.
 */

#include "multiPhysics/boussinesqThermalProcessor3D.h"
#include "multiPhysics/advectionDiffusion3D.h"
#include "multiPhysics/interparticlePotential.h"
#include "multiPhysics/shanChenLattices3D.h"
#include "multiPhysics/shanChenProcessor3D.h"
#include "multiPhysics/thermalDataAnalysis3D.h"
#include "multiPhysics/heLeeProcessor3D.h"
#include "multiPhysics/heLeeProcessor3D.h"
#include "multiPhysics/freeSurfaceUtil3D.h"
#include "multiPhysics/freeSurfaceModel3D.h"
#include "multiPhysics/freeSurfaceBoundaryCondition3D.h"
#include "multiPhysics/freeSurfaceInitializer3D.h"
#include "multiPhysics/freeSurfaceAnalysis3D.h"
#include "multiPhysics/multiFreeSurfaceModel3D.h"
#include "multiPhysics/createBubbles3D.h"
#include "multiPhysics/bubbleHistory3D.h"
#include "multiPhysics/bubbleMatch3D.h"
#include "multiPhysics/twoPhaseModel3D.h"
#include "multiPhysics/bodyForce3D.h"

