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

#include "palabos3D.h"
#include "palabos3D.hh"

#include <iostream>

using namespace plb;
typedef double T;

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    
    Array<Array<T,3>,3> M;
    M[0][0] = 3.; M[0][1] = 0.; M[0][2] = 0.;
    M[1][0] = 0.; M[1][1] = 7.; M[1][2] = 0.;
    M[2][0] = 0.; M[2][1] = 0.; M[2][2] = 8.;
    Array<Array<T,3>,3> V;
    Array<T,3> d;
    eigenDecomposition(M, V, d);
    std::cout << "Eigenvalue 0: " << d[0] << std::endl;
    std::cout << "Eigenvalue 1: " << d[1] << std::endl;
    std::cout << "Eigenvalue 2: " << d[2] << std::endl;

    std::cout << "Eigenvector 0: " << V[0][0] << ", " << V[1][0] << ", " << V[2][0] << std::endl;
    std::cout << "Eigenvector 1: " << V[0][1] << ", " << V[1][1] << ", " << V[2][1] << std::endl;
    std::cout << "Eigenvector 2: " << V[0][2] << ", " << V[1][2] << ", " << V[2][2] << std::endl;

    return 0;
}

