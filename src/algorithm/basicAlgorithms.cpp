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

#include "algorithm/basicAlgorithms.h"

namespace plb {

namespace algorithm {

std::vector<plint> primeFactor(plint value) {
    std::vector<plint> primeFactors;
    plint testFactor = 2;
    while (testFactor <= value) {
        if (value%testFactor==0) {
            value /= testFactor;
            primeFactors.push_back(testFactor);
        }
        else {
            ++testFactor;
        }
    }
    return primeFactors;
}

std::vector<plint> evenRepartition(plint value, plint d) {
    std::vector<plint> primeFactors = primeFactor(value);
    std::vector<plint> repartition(d);
    for (plint iRep=0; iRep<d; ++iRep) {
        repartition[iRep] = 1;
    }
    plint iDim=0;
    for (plint iPrime=(int)(primeFactors.size()-1); iPrime>=0; --iPrime) {
        repartition[iDim] *= primeFactors[iPrime];
        iDim = (iDim+1)%d;
    }
    return repartition;
}

}

}
