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

#include "gridRefinementUtil3D.hh"
#include "couplingInterfaceGenerator3D.hh"
#include "couplingActionsGenerator3D.hh"
#include "octree.hh"
#include "octreeGridGenerator.hh"
#include "refinementCriteria3D.hh"
#include "rescaleEngine.hh"
#include "boxLogic3D.hh"
#include "multiLevelScalarField3D.hh"
#include "multiLevelTensorField3D.hh"
#include "multiLevelNTensorField3D.hh"
#include "multiLevelWrapper3D.hh"
#include "multiLevelFieldGenerator3D.hh"
#include "gridRefinementFunctional3D.hh"
#include "dataAnalysisWrapper3D.hh"
