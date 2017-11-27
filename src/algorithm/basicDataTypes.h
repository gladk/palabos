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

#ifndef BASIC_DATA_TYPES_H
#define BASIC_DATA_TYPES_H

#include "core/globalDefs.h"
#include <set>
#include <typeinfo>

namespace plb {


template<typename T>
class EnumeratedUniqueObjects {
public:
    pluint operator[] (T const& object) {
        typename std::set<T>::iterator pos = objects.insert(object).first;
        return distance(objects.begin(), pos);
    }
private:
    std::set<T> objects;
};

struct TypeInfoComparator
{
    bool operator () (std::type_info const* t1, std::type_info const* t2) const {
        return t1->before(*t2);
    }
};


class EnumeratedUniqueTypeIds {
public:
    pluint operator[] (std::type_info const& object) {
        std::set<std::type_info const*, TypeInfoComparator>::iterator pos = objects.insert(&object).first;
        return std::distance(objects.begin(), pos);
    }
private:
    std::set<std::type_info const*, TypeInfoComparator> objects;
};

} // namespace plb

#endif  // BASIC_DATA_TYPES_H
