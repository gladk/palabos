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

#include "core/indexUtil.h"

namespace plb {

IndexCollection::IndexCollection(std::vector<plint> const& indexes_)
    : indexes(indexes_)
{ }

std::vector<plint> const& IndexCollection::get() const {
    return indexes;
}

IndexCollection operator&&(IndexCollection const& coll1, IndexCollection const& coll2) {
    std::vector<plint> result;
    std::vector<plint> const& ind1 = coll1.get();
    std::vector<plint> const& ind2 = coll2.get();

    std::vector<plint>::const_iterator iter_ind1 = ind1.begin();
    std::vector<plint>::const_iterator iter_ind2 = ind2.begin();

    while(!(iter_ind1==ind1.end() || iter_ind2==ind2.end())) {
        if (*iter_ind1 == *iter_ind2) {
            result.push_back(*iter_ind1);
            ++iter_ind1;
            ++iter_ind2;
        }
        else {
            if (*iter_ind1 < *iter_ind2) {
                ++iter_ind1;
            }
            else {
                ++iter_ind2;
            }
        }
    }
    return IndexCollection(result);
}

IndexCollection operator||(IndexCollection const& coll1, IndexCollection const& coll2)
{
    std::vector<plint> result;
    std::vector<plint> const& ind1 = coll1.get();
    std::vector<plint> const& ind2 = coll2.get();

    std::vector<plint>::const_iterator iter_ind1 = ind1.begin();
    std::vector<plint>::const_iterator iter_ind2 = ind2.begin();

    while(!(iter_ind1==ind1.end() && iter_ind2==ind2.end())) {
        if (iter_ind1==ind1.end()) {
            result.push_back(*iter_ind2);
            ++iter_ind2;
        }
        else if (iter_ind2==ind2.end()) {
            result.push_back(*iter_ind1);
            ++iter_ind1;
        }
        else {
            if (*iter_ind1==*iter_ind2) {
                result.push_back(*iter_ind1);
                ++iter_ind1;
                ++iter_ind2;
            }
            else {
                if (*iter_ind1<*iter_ind2) {
                    result.push_back(*iter_ind1);
                    ++iter_ind1;
                }
                else {
                    result.push_back(*iter_ind2);
                    ++iter_ind2;
                }
            }
        }
    }
    return IndexCollection(result);
}

IndexCollection operator!(IndexCollection const& coll) {
    std::vector<plint> result;
    std::vector<plint> const& ind = coll.get();

    std::vector<plint>::const_iterator iter_ind = ind.begin();

    plint negative_ind = 0;
    for (; iter_ind != ind.end(); ++iter_ind) {
        while (negative_ind < *iter_ind) {
            result.push_back(negative_ind++);
        }
        ++negative_ind;
    }
    return IndexCollection(result);
}
std::vector<plint> findIndexes(IndexCollection const& collection) {
    return collection.get();
}

}  // namespace plb

