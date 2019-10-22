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
 * Groups all the 3D include files in the complexDynamics directory.
*/

#include "complexDynamics/advectionDiffusionBoundaries.h"
#include "complexDynamics/advectionDiffusionBoundaryCondition3D.h"
#include "complexDynamics/adiabaticBoundaryProcessor3D.h"
#include "complexDynamics/advectionDiffusionDynamics.h"
#include "complexDynamics/advectionDiffusionProcessor3D.h"
#include "complexDynamics/advectionDiffusionUnits.h"
#include "complexDynamics/entropicDynamics.h"
#include "complexDynamics/mrtDynamics.h"
#include "complexDynamics/trtDynamics.h"
#include "complexDynamics/variableOmegaDynamics.h"
#include "complexDynamics/smagorinskyDynamics.h"
#include "complexDynamics/smagorinskyDynamics3D.h"
#include "complexDynamics/carreauDynamics.h"
#include "complexDynamics/carreauDynamicsTemplates.h"
#include "complexDynamics/carreauUnits.h"
#include "complexDynamics/carreauGlobalDefs.h"
#include "complexDynamics/dynamicSmagorinskyLattices3D.h"
#include "complexDynamics/dynamicSmagorinskyDynamics.h"
#include "complexDynamics/dynamicSmagorinskyProcessor3D.h"
#include "complexDynamics/externalForceMrtDynamics.h"
#include "complexDynamics/asinariModel.h"
#include "complexDynamics/wavePropagation.h"

