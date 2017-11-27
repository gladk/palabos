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
 * Groups all the include files for basic 3D dynamics.
 */
#include "core/globalDefs.h"
#include "core/plbComplex.h"
#include "core/plbInit.h"
#include "core/geometry3D.h"
#include "core/blockIdentifiers.h"
#include "core/dynamicsIdentifiers.h"
#include "core/units.h"
#include "core/dynamics.h"
#include "core/cell.h"
#include "core/blockStatistics.h"
#include "core/dataFieldBase3D.h"
#include "core/serializer.h"
#include "core/block3D.h"
#include "core/blockLatticeBase3D.h"
#include "core/latticeStatistics.h"
#include "core/plbTimer.h"
#include "core/plbRandom.h"
#include "core/plbLogFiles.h"
#include "core/indexUtil.h"
#include "core/multiBlockIdentifiers3D.h"
#include "core/nonLocalDynamics3D.h"
#include "core/functions.h"
#include "core/vectorFunction3D.h"
#include "core/realFunction3D.h"

