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
 * A timer class for generating random number -- header file.
 */
#ifndef PLB_RANDOM_H
#define PLB_RANDOM_H

#include <cstdlib>

namespace plb {

namespace global {

template<typename T>
struct PlbRandom {
    virtual ~PlbRandom() { }
    virtual void seed(pluint seedValue_) =0;
    virtual void iterate() =0;
    virtual T get(pluint position_) =0;
};

template<typename T>
inline T adjustRandomNumber(pluint rNumber) {
    return (T)rNumber;
}

template<>
inline double adjustRandomNumber<double>(pluint rNumber) {
    return (double)rNumber/(double)((pluint)RAND_MAX+1);
}

template<>
inline float adjustRandomNumber<float>(pluint rNumber) {
    return (float)rNumber/(float)((pluint)RAND_MAX+1);
}

template<typename T>
class DefaultRandom : public PlbRandom<T>
{
public:
    virtual void seed(pluint seedValue_) {
        seedValue = seedValue_;
    }
    virtual void iterate() {
        srand((unsigned int) seedValue);
        seedValue = (pluint)rand();
    }
    virtual T get(pluint position) {
        pluint currentSeed = seedValue+position;
        srand((unsigned int) (currentSeed%(pluint)RAND_MAX));
        for (pluint i=0; i<=(15+currentSeed/(pluint)RAND_MAX); ++i) {
            rand();
        }
        return adjustRandomNumber<T>((pluint)rand());
    }
private:
    pluint seedValue;
};

template<typename T>
inline PlbRandom<T>& plbRandom() {
    static DefaultRandom<T> defaultRandom;
    return defaultRandom;
}

}  // namespace global

}  // namespace plb

#endif  // PLB_RANDOM_H

