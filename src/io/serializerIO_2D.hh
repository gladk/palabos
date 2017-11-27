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

#ifndef SERIALIZER_IO_2D_HH
#define SERIALIZER_IO_2D_HH

#include "io/serializerIO_2D.h"
#include "core/plbDebug.h"
#include "core/runTimeDiagnostics.h"
#include "parallelism/mpiManager.h"
#include "multiBlock/multiBlockSerializer2D.h"
#include <istream>
#include <ostream>
#include <fstream>

namespace plb {

template<typename T>
std::ostream& block2ostream(std::ostream& ostr, Block2D const& block) {
    plint numDigits = 0; // Number of digits is handled by iostream manipulators,
                        // as opposed to being imposed here.
    serializerToAsciiStream<T>( block.getBlockSerializer (
                                    block.getBoundingBox(),
                                    global::IOpolicy().getIndexOrderingForStreams() ),
                                &ostr, numDigits );
    return ostr;
}

template<typename T>
std::istream& istream2block(std::istream& istr, Block2D& block) {
    asciiStreamToUnSerializer<T>( &istr,
                                  block.getBlockUnSerializer (
                                      block.getBoundingBox(),
                                          global::IOpolicy().getIndexOrderingForStreams()
                                  ) );
    return istr;
}

} // namespace plb

#endif  // SERIALIZER_IO_2D_HH
