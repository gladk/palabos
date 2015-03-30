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

#ifndef SPLINE_HH
#define SPLINE_HH

#include "core/runTimeDiagnostics.h"
#include "core/util.h"
#include "algorithm/spline.h"

#include <cstdio>
#include <limits>

namespace plb {

/* ***************** class Spline ***************************************** */

template<typename T>
Spline<T>::Spline(std::string fname)
{ 
    FILE *fp = fopen(fname.c_str(), "r");
    PLB_ASSERT(fp != NULL);

    long double tmp0, tmp1;
    while (!feof(fp))
        if (fscanf(fp, "%Lf%Lf", &tmp0, &tmp1) != EOF) {
            x.push_back((T) tmp0);
            y.push_back((T) tmp1);
        }

    fclose(fp);

    PLB_ASSERT(x.size() == y.size());
#ifdef PLB_DEBUG
    T eps = 100.0 * std::numeric_limits<T>::epsilon();
    for (plint i = 1; i < (plint) x.size(); i++)
        if (x[i] < x[i - 1] || util::fpequal(x[i], x[i - 1], eps)) {
                plbLogicError("The abscissa of the interpolation points must be monotonically increasing.");
        }
#endif // PLB_DEBUG
}

template<typename T>
Spline<T>::Spline(std::vector<T> const& x_, std::vector<T> const& y_)
    : x(x_),
      y(y_)
{ 
    PLB_ASSERT(x.size() == y.size());
#ifdef PLB_DEBUG
    T eps = 100.0 * std::numeric_limits<T>::epsilon();
    for (plint i = 1; i < (plint) x.size(); i++)
        if (x[i] < x[i - 1] || util::fpequal(x[i], x[i - 1], eps)) {
                plbLogicError("The abscissa of the interpolation points must be monotonically increasing.");
        }
#endif // PLB_DEBUG
}

/* ***************** class NaturalCubicSpline ***************************************** */

template<typename T>
NaturalCubicSpline<T>::NaturalCubicSpline(std::string fname)
    : CubicSpline<T>(fname)
{ 
    icache = 0;
    constructSpline();
}

template<typename T>
NaturalCubicSpline<T>::NaturalCubicSpline(std::vector<T> const& x_, std::vector<T> const& y_)
    : CubicSpline<T>(x_, y_)
{ 
    icache = 0;
    constructSpline();
}

template<typename T>
NaturalCubicSpline<T>* NaturalCubicSpline<T>::clone() const
{
    return new NaturalCubicSpline<T>(*this);
}

template<typename T>
void NaturalCubicSpline<T>::constructSpline()
{
    std::vector<T> const& x = this->getAbscissae();
    std::vector<T> const& y = this->getOrdinates();
    plint n = (plint) x.size();
    PLB_ASSERT( n>=2 );
    std::vector<T> a1(n-1), a2(n-1), a3(n), a4(n), a5(n);

    y1.resize(n-1);
    y2.resize(n);
    y3.resize(n-1);

    // Calculate the spline coefficients for polynomials of the form:
    // S_j(x) = y_j + y1_j * (x - x_j) + y2_j * (x - x_j)^2 + y3_j * (x - x_j)^3 

    for (plint i = 0; i < n - 1; i++)
        a1[i] = x[i + 1] - x[i];

    for (plint i = 1; i < n - 1; i++)
        a2[i] = 3.0 * (y[i + 1] - y[i]) / a1[i] - 3.0 * (y[i] - y[i - 1]) / a1[i - 1];

    a3[0] = 1.0;
    a4[0] = 0.0;
    a5[0] = 0.0;

    for (plint i = 1; i < n - 1; i++) {
        a3[i] = 2.0 * (x[i + 1] - x[i - 1]) - a1[i - 1] * a4[i - 1];
        a4[i] = a1[i] / a3[i];
        a5[i] = (a2[i] - a1[i - 1] * a5[i - 1]) / a3[i];
    }

    a3[n - 1] = 1.0;
    a5[n - 1] = 0.0;
    y2[n - 1] = 0.0;

    for (plint i = n - 2; i > -1; i--) {
        y2[i] = a5[i] - a4[i] * y2[i + 1];
        y1[i] = (y[i + 1] - y[i]) / a1[i] - a1[i] * (y2[i + 1] + 2.0 * y2[i]) / 3.0;
        y3[i] = (y2[i + 1] - y2[i]) / (3.0 * a1[i]);
    }

    y2.resize(n-1);
}

// return an index i in the range from il to ih for which x[i] <= t < x[i + 1]
template<typename T>
plint NaturalCubicSpline<T>::bsrch(T t, plint il, plint ih) const
{
    std::vector<T> const& x = this->getAbscissae();
    plint ilo = il;
    plint ihi = ih;

    plint i;
    while (ihi > (ilo + 1)) {
        i = (ihi + ilo) / 2;

        if (x[i] > t)
            ihi = i;
        else
            ilo = i;
    }

    return ilo;
}

template<typename T>
T NaturalCubicSpline<T>::getFunctionValue(T t) const
{
    std::vector<T> const& x = this->getAbscissae();
    std::vector<T> const& y = this->getOrdinates();
    plint n = (plint) x.size();
    PLB_ASSERT( n>=1 );

#ifdef PLB_DEBUG
    T x1 = x[0];
    T x2 = x[n-1];
    static T eps = 100.0 * std::numeric_limits<T>::epsilon();
    if ((t < x1 && !util::fpequal(t, x1, eps)) ||
        (t > x2 && !util::fpequal(t, x2, eps)))
        plbLogicError("Argument out of bounds.");
#endif // PLB_DEBUG

    if (t < x[icache] || t >= x[icache + 1]) {
        if (t > x[icache]) {
            icache = bsrch(t, icache, n - 1);
        } else {
            icache = bsrch(t, 0, icache);
        }
    }

    T xtmp = t - x[icache];
    T ytmp = y[icache] + xtmp * (y1[icache] + xtmp * (y2[icache] + xtmp * (y3[icache])));
    return ytmp;
}

template<typename T>
T NaturalCubicSpline<T>::getDerivativeValue(T t) const
{
    std::vector<T> const& x = this->getAbscissae();
    plint n = (plint) x.size();
    PLB_ASSERT( n>=1 );

#ifdef PLB_DEBUG
    T x1 = x[0];
    T x2 = x[n-1];

    static T eps = 100.0 * std::numeric_limits<T>::epsilon();
    if ((t < x1 && !util::fpequal(t, x1, eps)) ||
        (t > x2 && !util::fpequal(t, x2, eps)))
        plbLogicError("Argument out of bounds.");
#endif // PLB_DEBUG

    if (t < x[icache] || t >= x[icache + 1]) {
        if (t > x[icache]) {
            icache = bsrch(t, icache, n - 1);
        } else {
            icache = bsrch(t, 0, icache);
        }
    }

    T xtmp = t - x[icache];
    T ydtmp = y1[icache] + xtmp * (2.0 * y2[icache] + xtmp * (3.0 * y3[icache]));
    return ydtmp;
}

template<typename T>
T NaturalCubicSpline<T>::getSecondDerivativeValue(T t) const
{
    std::vector<T> const& x = this->getAbscissae();
    plint n = (plint) x.size();
    PLB_ASSERT( n>=1 );

#ifdef PLB_DEBUG
    T x1 = x[0];
    T x2 = x[n-1];

    static T eps = 100.0 * std::numeric_limits<T>::epsilon();
    if ((t < x1 && !util::fpequal(t, x1, eps)) ||
        (t > x2 && !util::fpequal(t, x2, eps)))
        plbLogicError("Argument out of bounds.");
#endif // PLB_DEBUG

    if (t < x[icache] || t >= x[icache + 1]) {
        if (t > x[icache]) {
            icache = bsrch(t, icache, n - 1);
        } else {
            icache = bsrch(t, 0, icache);
        }
    }

    T xtmp = t - x[icache];
    T yd2tmp = 2.0 * y2[icache] + 6.0 * y3[icache] * xtmp;
    return yd2tmp;
}

template<typename T>
T NaturalCubicSpline<T>::getThirdDerivativeValue(T t) const
{
    std::vector<T> const& x = this->getAbscissae();
    plint n = (plint) x.size();
    PLB_ASSERT( n>=1 );

#ifdef PLB_DEBUG
    T x1 = x[0];
    T x2 = x[n-1];

    static T eps = 100.0 * std::numeric_limits<T>::epsilon();
    if ((t < x1 && !util::fpequal(t, x1, eps)) ||
        (t > x2 && !util::fpequal(t, x2, eps)))
        plbLogicError("Argument out of bounds.");
#endif // PLB_DEBUG

    if (t < x[icache] || t >= x[icache + 1]) {
        if (t > x[icache]) {
            icache = bsrch(t, icache, n - 1);
        } else {
            icache = bsrch(t, 0, icache);
        }
    }

    T yd3tmp = 6.0 * y3[icache];
    return yd3tmp;
}

template<typename T>
T NaturalCubicSpline<T>::getIntegralValue() const
{
    std::vector<T> const& x = this->getAbscissae();
    std::vector<T> const& y = this->getOrdinates();
    plint n = (plint) x.size();

    T integral = (T) 0;
    for (plint i = 0; i < n-1; i++) {
        T dx = x[i + 1] - x[i];
        integral += dx * (y[i] + dx * (y1[i] / 2.0 + dx * (y2[i] / 3.0 + dx * (y3[i] / 4.0))));
    } 

    return integral;
}

template<typename T>
T NaturalCubicSpline<T>::getIntegralValue(T tmin, T tmax) const
{
    static T eps = 100.0 * std::numeric_limits<T>::epsilon();
    std::vector<T> const& x = this->getAbscissae();
    std::vector<T> const& y = this->getOrdinates();
    plint n = (plint) x.size();

#ifdef PLB_DEBUG
    T x1 = x[0];
    T x2 = x[n-1];

    if (tmin > tmax)
        plbLogicError("Invalid arguments.");
    if ((tmin < x1 && !util::fpequal(tmin, x1, eps)) ||
        (tmin > x2 && !util::fpequal(tmin, x2, eps)))
        plbLogicError("Argument out of bounds.");
    if ((tmax < x1 && !util::fpequal(tmax, x1, eps)) ||
        (tmax > x2 && !util::fpequal(tmax, x2, eps)))
        plbLogicError("Argument out of bounds.");
#endif // PLB_DEBUG

    if (util::fpequal(tmin, tmax, eps))
        return ((T) 0);

    if (tmin < x[icache] || tmin >= x[icache + 1]) {
        if (tmin > x[icache]) {
            icache = bsrch(tmin, icache, n - 1);
        } else {
            icache = bsrch(tmin, 0, icache);
        }
    }
    plint imin = icache;

    if (tmax < x[icache] || tmax >= x[icache + 1]) {
        if (tmax > x[icache]) {
            icache = bsrch(tmax, icache, n - 1);
        } else {
            icache = bsrch(tmax, 0, icache);
        }
    }
    plint imax = icache;

    plint i;
    T dxmin, dxmax, integral;

    if (imin == imax) {
        i = imin;
        dxmin = tmin - x[i];
        dxmax = tmax - x[i];
        integral = y[i] * (tmax - tmin) +
            (y1[i]/2.0) * (dxmax*dxmax - dxmin*dxmin) +
            (y2[i]/3.0) * (dxmax*dxmax*dxmax - dxmin*dxmin*dxmin) +
            (y3[i]/4.0) * (dxmax*dxmax*dxmax*dxmax - dxmin*dxmin*dxmin*dxmin);
        return integral;
    }

    i = imin;
    dxmin = tmin - x[i];
    dxmax = x[i + 1] - x[i];
    integral = y[i] * (x[i + 1] - tmin) +
        (y1[i]/2.0) * (dxmax*dxmax - dxmin*dxmin) +
        (y2[i]/3.0) * (dxmax*dxmax*dxmax - dxmin*dxmin*dxmin) +
        (y3[i]/4.0) * (dxmax*dxmax*dxmax*dxmax - dxmin*dxmin*dxmin*dxmin);

    T dx;
    for (i = imin + 1; i < imax; i++) {
        dx = x[i + 1] - x[i];
        integral += dx * (y[i] + dx * (y1[i] / 2.0 + dx * (y2[i] / 3.0 + dx * (y3[i] / 4.0))));
    } 

    i = imax;
    dx = tmax - x[i];
    integral += dx * (y[i] + dx * (y1[i] / 2.0 + dx * (y2[i] / 3.0 + dx * (y3[i] / 4.0))));

    return integral;
}

}  // namespace plb

#endif  // SPLINE_HH
