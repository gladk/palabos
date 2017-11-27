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

#include "io/base64.h"
#include "io/base64.hh"
#include "io/serializerIO_2D.h"
#include "io/serializerIO_2D.hh"
#include "io/serializerIO_3D.h"
#include "io/serializerIO_3D.hh"

namespace plb {

// All of the following is a workaround to the following problem: on a
// 32-bit machine where pluint plint is the same as pluint, Base64Encoder
// needs to be instantiated on one integer type only. On some 64-bit
// platforms however, pluint is not equal to pluint int. In that case,
// Base64Encoder needs to be instantiated twice. It is however not
// possible to instantiate this class first on pluint plint and then
// on pluint, because this yields a double instantiation, and thus
// an error, where these types are the same. To avoid this problem,
// the chosen instantiation types are pluint plint and pluint where
// these types are different, and truc and pluint plint where
// they are similar. A template-based if-then-else construct is used
// to distinguish the two cases.

template<bool areEqual> struct DistinctUint;

template<>
struct DistinctUint<true> {
    typedef unsigned char T1;
    typedef pluint T2;
};

template<>
struct DistinctUint<false> {
    typedef unsigned int T1;
    typedef pluint T2;
};

typedef DistinctUint<sizeof(unsigned int)==sizeof(pluint)>::T1 T1;
typedef DistinctUint<sizeof(unsigned int)==sizeof(pluint)>::T2 T2;

template class Base64Encoder<FLOAT_T>;
template class Base64Encoder<T1>;
template class Base64Encoder<T2>;

template class Base64Decoder<FLOAT_T>;
template class Base64Decoder<T1>;
template class Base64Decoder<T2>;

}
