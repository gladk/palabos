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

#ifndef PLB_DEBUG_H
#define PLB_DEBUG_H

#include <cassert>

#ifdef PLB_DEBUG

    #define PLB_ASSERT( COND )        assert( COND );
    #define PLB_PRECONDITION( COND )  assert( COND );
    #define PLB_POSTCONDITION( COND ) assert( COND );
    #define PLB_STATECHECK( A,B )     assert( (A) == (B) );

#else

    #define PLB_ASSERT( COND )
    #define PLB_PRECONDITION( COND )
    #define PLB_POSTCONDITION( COND )
    #define PLB_STATECHECK( A,B )

#endif  // PLB_DEBUG

namespace plb {

// Programmatically enable core dumps for POSIX systems.
// In a parallel program, this function must be called after plbInit.
void enableCoreDumps();

// Make stdout and stderr unbuffered for better debugging.
void unbufferOutputStdStreams();

}  // namespace plb

#endif  // PLB_DEBUG_H
