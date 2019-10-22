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
 * Geometric operations on collections of 2D domains -- header file.
 */
#ifndef DOMAIN_MANIPULATION_2D_H
#define DOMAIN_MANIPULATION_2D_H

#include "core/geometry2D.h"
#include <vector>

namespace plb {

struct DomainAndId2D {
    DomainAndId2D(Box2D domain_, plint id_)
        : domain(domain_),
          id(id_)
    { }
    Box2D domain;
    plint id;
};

/// Compute mutual intersections between domains, and remove overlaps.
/** The union of all domains in the resulting vector is the same as the union
 *  of all domains in the argument, but the resulting domains never overlap.
 */
std::vector<DomainAndId2D> getNonOverlapingBlocks(std::vector<DomainAndId2D> const& domainsWithId);

/// Compute common intersections among several collections of domains.
/** For each intersection, the IDs of the original domain in the original
 *  collection is also stored.
 */
void intersectDomainsAndIds(std::vector<std::vector<DomainAndId2D> > const& domainsWithId,
                            std::vector<Box2D>& finalDomains,
                            std::vector<std::vector<plint> >& finalIds);

} // namespace plb

#endif  // DOMAIN_MANIPULATION_2D_H
