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
 * LB initialisation routine -- implementation file.
 */

#include "core/plbInit.h"
#include "core/plbInit.hh"
#include "core/plbProfiler.h"
#include "core/plbRandom.h"
#include "core/runTimeDiagnostics.h"
#include <sstream>

namespace plb {

void plbInit(int *argc, char ***argv, bool verbous) {
    global::mpi().init(argc, argv, verbous);
    global::mainArguments().setArgs(*argc, argv);
    global::plbRandom<float>().seed(10);
    global::plbRandom<double>().seed(10);
    global::plbRandom<plint>().seed(10);
}

void plbInit() {
    global::mpi().init();
    global::plbRandom<float>().seed(10);
    global::plbRandom<double>().seed(10);
    global::plbRandom<plint>().seed(10);
}

namespace global {

    MainArgv::MainArgv(std::string argument_, int whichArg_)
        : argument(argument_),
          whichArg(whichArg_)
    { }

    MainArgv:: operator std::string() const {
        return argument;
    }

    int MainArgs::argc() const {
        if (arguments.empty()) {
            plbLogicError("Can't access command-line arguments: they have not yet been initialized.");
        }
        return (int) arguments.size();
    }

    MainArgv MainArgs::argv(int whichArg) const {
        if (arguments.empty()) {
            plbLogicError("Can't access command-line arguments: they have not yet been initialized.");
        }
        if (whichArg >= argc()) {
            std::stringstream errMessage;
            errMessage << "Can\'t read command-line argument " << whichArg << ": ";
            if (argc()==1) {
                errMessage << "there are no command-line arguments";
            }
            else if (argc()==2) {
                errMessage << "there is only one command-line argument";
            }
            else {
                errMessage << "there are only " << argc()-1 << " arguments";
            }
            plbIOError(errMessage.str());
        }
        return MainArgv(arguments[whichArg], whichArg);
    }

    MainArgs::MainArgs()
    { }

    void MainArgs::setArgs(int argcValue, char*** argvPointer) {
        arguments.resize(argcValue);
        for (plint iArg=0; iArg<argcValue; ++iArg) {
            arguments[iArg] = std::string((*argvPointer)[iArg]);
        }
    }

    void MainArgs::setArgs(std::vector<std::string> arguments_) {
        arguments = arguments_;
    }

    MainArgs& mainArguments() {
        static MainArgs instance;
        return instance;
    }

    int argc() {
        return mainArguments().argc();
    }

    MainArgv argv(int whichArg) {
        return mainArguments().argv(whichArg);
    }

    template void MainArgv::read<int>(int& variable);
    template void MainArgv::read<double>(double& variable);
    template void MainArgv::read<std::string>(std::string& variable);

}  // namespace global

}  // namespace plb
