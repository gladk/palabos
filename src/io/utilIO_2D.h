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
 * I/O utilities -- header file.
 */

#ifndef UTIL_IO_2D_H
#define UTIL_IO_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/atomicContainerBlock2D.h"
#include "multiBlock/multiBlock2D.h"
#include "io/plbFiles.h"

namespace plb {

/* Save the current state of the simulation for restarting. */
void saveState(std::vector<MultiBlock2D*> blocks, plint iteration, bool saveDynamicContent,
        FileName xmlFileName, FileName baseFileName, plint fileNamePadding = 8);

/* Load the state of the simulation from checkpoint files for restarting. */
void loadState(std::vector<MultiBlock2D*> blocks, plint& iteration, bool saveDynamicContent,
        FileName xmlFileName);

/* Check for user-driven execution abortion, and save the state of the simulation. */
bool abortExecution(FileName abortFileName, std::vector<MultiBlock2D*> blocks, plint iteration,
        bool saveDynamicContent, FileName xmlFileName, FileName baseFileName,
        plint fileNamePadding = 8);

}  // namespace plb

#endif  // UTIL_IO_2D_H

