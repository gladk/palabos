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

#ifndef INTERPOLATION_2D_HH
#define INTERPOLATION_2D_HH

#include "core/util.h"
#include "core/globalDefs.h"
#include "algorithm/interpolation2D.h"

#include <cstdio>
#include <limits>

namespace plb {

/* ***************** class Interpolation2D ***************************************** */

template<typename T>
Interpolation2D<T>::Interpolation2D(std::string fname)
{ 
    FILE *fp = fopen(fname.c_str(), "r");
    PLB_ASSERT(fp != 0);

    long int tmpNx, tmpNy;
    plint nx, ny;
    if (fscanf(fp, "%ld%ld", &tmpNx, &tmpNy) == 2) {
        nx = (plint) tmpNx;
        ny = (plint) tmpNy;
    } else {
        PLB_ASSERT(false);
    }
    PLB_ASSERT(nx >= 2);
    PLB_ASSERT(ny >= 2);
    x.resize(nx);
    y.resize(ny);
    f.resize(nx);
    for (plint i = 0; i < nx; i++) {
        f[i].resize(ny);
    }

    long double tmp;
    for (plint i = 0; i < nx; i++) {
        if (fscanf(fp, "%Lf", &tmp) == 1) {
            x[i] = (T) tmp;
        } else {
            PLB_ASSERT(false);
        }
    }
    for (plint i = 0; i < ny; i++) {
        if (fscanf(fp, "%Lf", &tmp) == 1) {
            y[i] = (T) tmp;
        } else {
            PLB_ASSERT(false);
        }
    }
    for (plint i = 0; i < nx; i++) {
        for (plint j = 0; j < ny; j++) {
            if (fscanf(fp, "%Lf", &tmp) == 1) {
                f[i][j] = (T) tmp;
            } else {
                PLB_ASSERT(false);
            }
        }
    }

    fclose(fp);

#ifdef PLB_DEBUG
    for (plint i = 1; i < nx; i++) {
        PLB_ASSERT(util::greaterThan(x[i], x[i - 1]));
    }
    for (plint i = 1; i < ny; i++) {
        PLB_ASSERT(util::greaterThan(y[i], y[i - 1]));
    }
#endif
}

template<typename T>
Interpolation2D<T>::Interpolation2D(std::vector<T> const& x_, std::vector<T> const& y_, std::vector<std::vector<T> > const& f_)
    : x(x_),
      y(y_),
      f(f_)
{ 
    PLB_ASSERT(x.size() >= 2);
    PLB_ASSERT(y.size() >= 2);
    PLB_ASSERT(f.size() == x.size());
#ifdef PLB_DEBUG
    for (plint i = 0; i < (plint) x.size(); i++) {
        PLB_ASSERT(f[i].size() == y.size());
    }
    for (plint i = 1; i < (plint) x.size(); i++) {
        PLB_ASSERT(util::greaterThan(x[i], x[i - 1]));
    }
    for (plint i = 1; i < (plint) y.size(); i++) {
        PLB_ASSERT(util::greaterThan(y[i], y[i - 1]));
    }
#endif
}

template<typename T>
void Interpolation2D<T>::transformCoordinatesX(T scale, T offset)
{
    for (plint i = 0; i < (plint) x.size(); i++) {
        x[i] = scale * x[i] + offset;
    }
}

template<typename T>
void Interpolation2D<T>::transformCoordinatesY(T scale, T offset)
{
    for (plint i = 0; i < (plint) y.size(); i++) {
        y[i] = scale * y[i] + offset;
    }
}

template<typename T>
void Interpolation2D<T>::transformFunctionValues(T scale, T offset)
{
    for (plint i = 0; i < (plint) x.size(); i++) {
        for (plint j = 0; j < (plint) y.size(); j++) {
            f[i][j] = scale * f[i][j] + offset;
        }
    }
}

template<typename T>
void Interpolation2D<T>::exportInXYZFormat(std::string fname, int numDecimalDigits) const
{ 
    if (global::mpi().isMainProcessor()) {
        FILE *fp = fopen(fname.c_str(), "w");
        PLB_ASSERT(fp != 0);
        int d = numDecimalDigits;
        for (plint i = 0; i < (plint) x.size(); i++) {
            for (plint j = 0; j < (plint) y.size(); j++) {
                double xd = (double) x[i];
                double yd = (double) y[j];
                double fd = (double) f[i][j];
                fprintf(fp, "% .*e % .*e % .*e\n", d, xd, d, yd, d, fd);
            }
        }
        fclose(fp);
    }
}


/* ***************** class LinearInterpolation2D ***************************************** */

template<typename T>
LinearInterpolation2D<T>* LinearInterpolation2D<T>::clone() const
{
    return new LinearInterpolation2D<T>(*this);
}

// Return an index i in the range from il to ih for which v[i] <= t < v[i + 1]
template<typename T>
plint LinearInterpolation2D<T>::bsrch(std::vector<T> const& v, T t, plint il, plint ih) const
{
    plint ilo = il;
    plint ihi = ih;

    plint i;
    while (ihi > (ilo + 1)) {
        i = (ihi + ilo) / 2;

        if (v[i] > t)
            ihi = i;
        else
            ilo = i;
    }

    return ilo;
}

template<typename T>
void LinearInterpolation2D<T>::locate(std::vector<T> const& v, T& t, plint& cache) const
{
    plint n = (plint) v.size();
    if (util::lessEqual(t, v[0])) {
        t = v[0];
        cache = 0;
    } else if (util::greaterEqual(t, v[n - 1])) {
        t = v[n - 1];
        cache = n - 2;
    } else {
        if (t < v[cache]) {
            cache = bsrch(v, t, 0, cache);
        } else if (t >= v[cache + 1]) {
            cache = bsrch(v, t, cache, n - 1);
        }
    }
}

template<typename T>
T LinearInterpolation2D<T>::getFunctionValue(T x0, T y0) const
{
    std::vector<T> const& x = this->getCoordinatesX();
    std::vector<T> const& y = this->getCoordinatesY();
    std::vector<std::vector<T> > const& f = this->getFunctionValues();

    locate(x, x0, icache);
    locate(y, y0, jcache);

    T x1 = x[icache    ];
    T x2 = x[icache + 1];
    T y1 = y[jcache    ];
    T y2 = y[jcache + 1];

    T f11 = f[icache    ][jcache    ];
    T f12 = f[icache    ][jcache + 1];
    T f21 = f[icache + 1][jcache    ];
    T f22 = f[icache + 1][jcache + 1];

    T ftmp = (f11 * (x2 - x0) * (y2 - y0) +
              f12 * (x2 - x0) * (y0 - y1) +
              f21 * (x0 - x1) * (y2 - y0) +
              f22 * (x0 - x1) * (y0 - y1)) / ((x2 - x1) * (y2 - y1));
    return ftmp;
}

}  // namespace plb

#endif  // INTERPOLATION_2D_HH
