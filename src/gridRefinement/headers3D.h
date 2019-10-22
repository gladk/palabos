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

#include "gridRefinementUtil3D.h"
#include "couplingInterfaceGenerator3D.h"
#include "couplingActionsGenerator3D.h"
#include "octree.h"
#include "octreeGridStructure.h"
#include "octreeGridGenerator.h"
#include "refinementCriteria3D.h"
#include "rescaleEngine.h"
#include "boxLogic3D.h"
#include "multiLevel3D.h"
#include "multiLevelScalarField3D.h"
#include "multiLevelTensorField3D.h"
#include "multiLevelNTensorField3D.h"
#include "multiLevelWrapper3D.h"
#include "multiLevelFieldGenerator3D.h"
#include "gridRefinementFunctional3D.h"
#include "dataAnalysisWrapper3D.h"

