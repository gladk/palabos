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

#ifndef TIME_PERIODIC_SIGNAL_H
#define TIME_PERIODIC_SIGNAL_H

#include <string>
#include "algorithm/spline.h"

namespace plb {

template<typename T>
class TimePeriodicSignal {
public:
    TimePeriodicSignal() : period() { }
    TimePeriodicSignal(std::string fname);
    TimePeriodicSignal(std::vector<T> const& t, std::vector<T> const& x);
    ~TimePeriodicSignal() { }
    TimePeriodicSignal<T>* clone() const;
    std::vector<T> const& getTimeValues() const { return signal.getAbscissae(); }
    std::vector<T>& getTimeValues() { return signal.getAbscissae(); }
    std::vector<T> const& getSignalValues() const { return signal.getOrdinates(); }
    std::vector<T>& getSignalValues() { return signal.getOrdinates(); }
    T getSignalValue(T t) const;
    T getDerivativeValue(T t) const;
    T getSecondDerivativeValue(T t) const;
    T getThirdDerivativeValue(T t) const;
    T getIntegralValue() const;
    T getIntegralValue(T tmin, T tmax) const;
private:
    NaturalCubicSpline<T> signal;
    T period;
};

}  // namespace plb

#endif  // TIME_PERIODIC_SIGNAL_H

