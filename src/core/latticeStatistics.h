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


#ifndef LATTICE_STATISTICS_H
#define LATTICE_STATISTICS_H

#include "core/globalDefs.h"
#include "core/blockStatistics.h"

namespace plb {

namespace LatticeStatistics {
    enum { avRhoBar=0, avUSqr=1, maxUSqr=0 };
}

template<typename T>
inline void gatherStatistics(BlockStatistics& statistics, const T& rhoBar, const T& uSqr) {
    statistics.gatherAverage(LatticeStatistics::avRhoBar, static_cast<double>(rhoBar));
    statistics.gatherAverage(LatticeStatistics::avUSqr, static_cast<double>(uSqr));
    statistics.gatherMax(LatticeStatistics::maxUSqr, static_cast<double>(uSqr));
    statistics.incrementStats();
}

} // namespace plb

#endif
