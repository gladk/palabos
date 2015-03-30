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

#ifndef NONLINEAR_EQUATION_SOLVERS_HH
#define NONLINEAR_EQUATION_SOLVERS_HH

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "algorithm/nonlinearEquationSolvers.h"

#include <cmath>

namespace plb {

/* ***************** class NewtonRaphsonMethod ***************************************** */

template<typename T>
NewtonRaphsonMethod<T>::NewtonRaphsonMethod(T tolerance_, plint maxNumOfIterations_)
    : tolerance(tolerance_),
      maxNumOfIterations(maxNumOfIterations_)
{
    PLB_ASSERT(tolerance > 0.0);
    PLB_ASSERT(maxNumOfIterations > 0);
}

template<typename T>
NewtonRaphsonMethod<T>* NewtonRaphsonMethod<T>::clone() const
{
    return new NewtonRaphsonMethod<T>(*this);
}

template<typename T>
bool NewtonRaphsonMethod<T>::solve(Function f, Function df, void *data, T initialGuess, T& solution, T& absoluteError) const
{
    bool converged = false;
    T xOld = initialGuess;

    T xNew = xOld;
    T error = -1.0;
    for (plint i = 0; i < maxNumOfIterations; i++) {
        xNew = xOld - f(xOld, data) / df(xOld, data);
        error = std::fabs(xNew - xOld);
        if (error <= tolerance) {
            converged = true;
            break;
        }
        xOld = xNew;
    }

    solution = xNew;
    absoluteError = error;
    return(converged);
}

/* ***************** Bisection Method ***************************************** */

template<typename T, class Function>
bool bisect(Function const& function, T x0, T x1, T xacc, plint maxIter, T& result)
{
    PLB_ASSERT(maxIter > 0);
    PLB_ASSERT(xacc > T());

    T dx, f, fmid, xmid;

    f=function(x0);
    fmid=function(x1);
    if (f*fmid >= T()) {
        result = T();
        return false;
    }

    // Orient the search so that f>0 lies at x+dx.
    if (f<T()) {
        dx = x1-x0;
        result = x0;
    }
    else {
        dx = x0-x1;
        result = x1;
    }

    for (int j=0; j<maxIter; ++j) { // Bisection loop.
        dx *= 0.5;
        xmid = result + dx;
        fmid = function(xmid);

        if (fmid <= T()) {
            result=xmid;
        }

        if (std::fabs(dx) < xacc || fmid == T()) {
            return true;
        }
    }

    result = T();
    return false;
}

}  // namespace plb

#endif  // NONLINEAR_EQUATION_SOLVERS_HH
