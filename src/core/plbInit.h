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
 * LB initialisation routine -- header file.
 */
#ifndef PLB_INIT_H
#define PLB_INIT_H

#include "core/globalDefs.h"
#include "parallelism/mpiManager.h"
#include "io/parallelIO.h"
#include <vector>
#include <string>

namespace plb {

void plbInit(int *argc, char ***argv, bool verbous=false);
void plbInit();

namespace global {

    class MainArgv {
    public:
        MainArgv(std::string argument_, int whichArg_);
        operator std::string() const;
        template<typename T> void read(T& variable);
        template<typename T> bool readNoThrow(T& variable);
    private:
        std::string argument;
        int whichArg;
    };

    class MainArgs {
    public:
        int argc() const;
        MainArgv argv(int whichArg) const;
        void setArgs(int argcValue, char*** argvPointer);
        void setArgs(std::vector<std::string> arguments_);
    private:
        MainArgs();
    private:
        std::vector<std::string> arguments;

    friend MainArgs& mainArguments();
    };


    MainArgs& mainArguments();
    int argc();
    MainArgv argv(int whichArg);

}  // namespace global

}  // namespace plb

#endif  // PLB_INIT_H
