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
 * LB initialisation routine -- generic code.
 */
#ifndef PLB_INIT_HH
#define PLB_INIT_HH

#include "core/plbInit.h"
#include "core/runTimeDiagnostics.h"
#include <sstream>

namespace plb {

namespace global {

template<typename T>
void MainArgv::read(T& variable) {
    T tmp = T();
    std::stringstream argStream;
    argStream << argument;
    argStream >> tmp;
    if (!argStream) {
        std::stringstream message;
        message << "Problem reading command-line argument " << whichArg;
        plbIOError(message.str());
    }
    variable = tmp;
}

template<typename T>
bool MainArgv::readNoThrow(T& variable) {
    T tmp = T();
    std::stringstream argStream;
    argStream << argument;
    argStream >> tmp;
    if (!argStream) {
        return false;
    }
    variable = tmp;
    return true;
}

}  // namespace global

}  // namespace plb

#endif  // PLB_INIT_HH
