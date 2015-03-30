/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
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

#ifndef TIME_PERIODIC_SIGNAL_HH
#define TIME_PERIODIC_SIGNAL_HH

#include "core/runTimeDiagnostics.h"
#include "core/util.h"
#include "algorithm/spline.h"
#include "algorithm/timePeriodicSignal.h"

#include <limits>

namespace plb {

template<typename T>
TimePeriodicSignal<T>::TimePeriodicSignal(std::string fname)
    : signal(fname)
{
    std::vector<T> const& t = signal.getAbscissae();
    std::vector<T> const& x = signal.getOrdinates();
    plint n = (plint) t.size();

#ifdef PLB_DEBUG
    T eps = (T) 100.0 * std::numeric_limits<T>::epsilon();
    if (!util::fpequal(t[0], (T) 0, eps))
        plbLogicError("The first given value of the signal must be at time 0.0.");
    if (!util::fpequal(x[0], x[n-1], eps))
        plbLogicError("The last given value of the signal must be equal to the first one.");
#endif // PLB_DEBUG

    period = t[n-1] - t[0];

#ifdef PLB_DEBUG
    if (util::fpequal(period, (T) 0, eps))
        plbLogicError("The period must be greater than zero.");
#endif // PLB_DEBUG
}

template<typename T>
TimePeriodicSignal<T>::TimePeriodicSignal(std::vector<T> const& t, std::vector<T> const& x)
    : signal(t, x)
{
    plint n = (plint) t.size();
    PLB_ASSERT(n >= 1);

#ifdef PLB_DEBUG
    T eps = (T) 100.0 * std::numeric_limits<T>::epsilon();
    if (!util::fpequal(t[0], (T) 0, eps))
        plbLogicError("The first given value of the signal must be at time 0.0.");
    if (!util::fpequal(x[0], x[n-1], eps))
        plbLogicError("The last given value of the signal must be equal to the first one.");
#endif // PLB_DEBUG

    period = t[n-1] - t[0];

#ifdef PLB_DEBUG
    if (util::fpequal(period, (T) 0, eps))
        plbLogicError("The period must be greater than zero.");
#endif // PLB_DEBUG
}

template<typename T>
TimePeriodicSignal<T>* TimePeriodicSignal<T>::clone() const {
    return new TimePeriodicSignal<T>(*this);
}

template<typename T>
T TimePeriodicSignal<T>::getSignalValue(T t) const
{
#ifdef PLB_DEBUG
    static T eps = (T) 100.0 * std::numeric_limits<T>::epsilon();
    if (t < (T) 0 && !util::fpequal(t, (T) 0, eps))
        plbLogicError("The time value must be greater equal to 0.0.");
#endif // PLB_DEBUG

    T trel = t - (plint) (t/period) * period;
    return signal.getFunctionValue(trel);
}

template<typename T>
T TimePeriodicSignal<T>::getDerivativeValue(T t) const
{
#ifdef PLB_DEBUG
    static T eps = (T) 100.0 * std::numeric_limits<T>::epsilon();
    if (t < (T) 0 && !util::fpequal(t, (T) 0, eps))
        plbLogicError("The time value must be greater equal to 0.0.");
#endif // PLB_DEBUG

    T trel = t - (plint) (t/period) * period;
    return signal.getDerivativeValue(trel);
}

template<typename T>
T TimePeriodicSignal<T>::getSecondDerivativeValue(T t) const
{
#ifdef PLB_DEBUG
    static T eps = (T) 100.0 * std::numeric_limits<T>::epsilon();
    if (t < (T) 0 && !util::fpequal(t, (T) 0, eps))
        plbLogicError("The time value must be greater equal to 0.0.");
#endif // PLB_DEBUG

    T trel = t - (plint) (t/period) * period;
    return signal.getSecondDerivativeValue(trel);
}

template<typename T>
T TimePeriodicSignal<T>::getThirdDerivativeValue(T t) const
{
#ifdef PLB_DEBUG
    static T eps = (T) 100.0 * std::numeric_limits<T>::epsilon();
    if (t < (T) 0 && !util::fpequal(t, (T) 0, eps))
        plbLogicError("The time value must be greater equal to 0.0.");
#endif // PLB_DEBUG

    T trel = t - (plint) (t/period) * period;
    return signal.getThirdDerivativeValue(trel);
}

template<typename T>
T TimePeriodicSignal<T>::getIntegralValue() const
{
    return signal.getIntegralValue();
}

template<typename T>
T TimePeriodicSignal<T>::getIntegralValue(T tmin, T tmax) const
{
    static T eps = (T) 100.0 * std::numeric_limits<T>::epsilon();

#ifdef PLB_DEBUG
    if (tmin > tmax)
        plbLogicError("Invalid arguments.");
    if (tmin < (T) 0 && !util::fpequal(tmin, (T) 0, eps))
        plbLogicError("Invalid arguments.");
    if (tmax < (T) 0 && !util::fpequal(tmax, (T) 0, eps))
        plbLogicError("Invalid arguments.");
#endif // PLB_DEBUG

    if (util::fpequal(tmin, tmax, eps))
        return ((T) 0);

    T dt = tmax - tmin;
    plint k = (plint) (dt / period);

    T integral = (k == 0) ? (T) 0 : k * signal.getIntegralValue();

    T tmin_rel = tmin - (plint) (tmin/period) * period;
    T tmax_rel = tmax - (plint) (tmax/period) * period;

    integral += signal.getIntegralValue(tmin_rel, tmax_rel);

    return integral;
}

}  // namespace plb

#endif  // TIME_PERIODIC_SIGNAL_HH
