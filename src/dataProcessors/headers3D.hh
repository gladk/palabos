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
 * Groups all the 3D generic implementation files for the directory "dataProcessors".
 */
#include "dataProcessors/dataAnalysisFunctional3D.hh"
#include "dataProcessors/dataAnalysisWrapper3D.hh"
#include "dataProcessors/dataInitializerFunctional3D.hh"
#include "dataProcessors/dataInitializerWrapper3D.hh"
#include "dataProcessors/metaStuffFunctional3D.hh"
#include "dataProcessors/metaStuffWrapper3D.hh"
#include "dataProcessors/ntensorAnalysisFunctional3D.hh"
#include "dataProcessors/ntensorAnalysisWrapper3D.hh"
// Include 2D versions, because they are required, for example to save 2D
// images from 3D data.
#include "dataProcessors/dataAnalysisFunctional2D.hh"
#include "dataProcessors/dataAnalysisWrapper2D.hh"
#include "dataProcessors/dataInitializerFunctional2D.hh"
#include "dataProcessors/dataInitializerWrapper2D.hh"
#include "dataProcessors/metaStuffWrapper2D.hh"

