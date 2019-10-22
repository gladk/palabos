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

#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "core/globalDefs.h"

#include <vector>

namespace plb {

template<typename T>
class GaussLegendreQuadrature {
public:
    typedef T (*IntegralKernel)(T x, void *data);
public:
    GaussLegendreQuadrature(plint n_, plint maxNumOfIterations_=64);
    GaussLegendreQuadrature<T>* clone() const;
    T evaluateIntegral(IntegralKernel integralKernel, void *data, T x0, T x1) const;
    plint getN() const { return n; }
    std::vector<T> const& getNodes() const { return nodes; };
    std::vector<T> const& getWeights() const { return weights; };
private:
    void evaluateLegendrePolynomialAndDerivative(plint k, T x, T& L_k, T& dL_k) const;
    void evaluateNodesAndWeights();
private:
    plint n;
    std::vector<T> nodes;
    std::vector<T> weights;
    plint maxNumOfIterations;
};

}  // namespace plb

#endif  // QUADRATURE_H

