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

#include "io/serializerIO_2D.h"
#include "core/plbDebug.h"
#include "core/runTimeDiagnostics.h"
#include "parallelism/mpiManager.h"
#include "multiBlock/multiBlockSerializer2D.h"
#include <istream>
#include <ostream>
#include <fstream>

namespace plb {

void saveBinaryBlock(Block2D const& block, std::string fName, bool enforceUint) {
    std::ofstream* ostr = 0;
    bool isOK = true;
    if (global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fName.c_str());
        isOK = (bool)(*ostr);
    }
    plbMainProcIOError( !isOK, std::string("Could not open binary file ")+
                               fName+std::string(" for saving") );
    serializerToBase64Stream( block.getBlockSerializer(block.getBoundingBox(),
                              global::IOpolicy().getIndexOrderingForStreams()), ostr, enforceUint );
    delete ostr;
}

void loadBinaryBlock(Block2D& block, std::string fName, bool enforceUint) {
    std::ifstream* istr = 0;
    bool isOK = true;
    if (global::mpi().isMainProcessor()) {
        istr = new std::ifstream(fName.c_str());
        isOK = (bool)(*istr);
    }
    plbMainProcIOError( !isOK, std::string("Could not open binary file ")+
                               fName+std::string(" for reading") );
    base64StreamToUnSerializer (
            istr, block.getBlockUnSerializer (
                block.getBoundingBox(),
                global::IOpolicy().getIndexOrderingForStreams()),
            enforceUint );
    delete istr;
}

} // namespace plb

