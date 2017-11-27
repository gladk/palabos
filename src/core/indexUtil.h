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
 * Templates for finding indexes for a specific subset of the neighborhood
 *  -- header file
 */
#ifndef INDEX_UTIL_H
#define INDEX_UTIL_H

#include "core/globalDefs.h"
#include <vector>

namespace plb {

class IndexCollection {
public:
    /// Assumption: indexes are sorted in asceding order.
    IndexCollection(std::vector<plint> const& indexes_);
    std::vector<plint> const& get() const;
private:
    std::vector<plint> indexes;
};

IndexCollection operator&&(IndexCollection const& coll1, IndexCollection const& coll2);

IndexCollection operator||(IndexCollection const& coll1, IndexCollection const& coll2);

IndexCollection operator!(IndexCollection const& coll);

template<typename T, template<typename U> class Descriptor>
class Index {
public:
    Index(plint direction_);
    IndexCollection operator==(plint value) const;
    IndexCollection operator<(plint value) const;
    IndexCollection operator<=(plint value) const;
    IndexCollection operator>(plint value) const;
    IndexCollection operator>=(plint value) const;
private:
    plint direction;
};

std::vector<plint> findIndexes(IndexCollection const& collection);

}  // namespace plb

#endif  // INDEX_UTIL_H
