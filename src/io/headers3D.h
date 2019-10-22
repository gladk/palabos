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
 * Groups all the 3D include files for the directory io.
 */

#include "io/base64.h"
#include "io/serializerIO.h"
#include "io/serializerIO_3D.h"
#include "io/vtkDataOutput.h"
#include "io/sparseVtkDataOutput.h"
#include "io/vtkStructuredDataOutput.h"
#include "io/parallelIO.h"
#include "io/colormaps.h"
#include "io/imageWriter.h"
#include "io/endianness.h"
#include "io/plbFiles.h"
#include "io/multiBlockReader3D.h"
#include "io/multiBlockWriter3D.h"
#include "io/utilIO_3D.h"
#include "io/transientStatistics3D.h"

