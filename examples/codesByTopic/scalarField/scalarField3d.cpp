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
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace std;
typedef double T;


int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    const int N = 300;
    ImageWriter<T> imageWriter("leeloo");

    MultiScalarField3D<T> xCoord(N+1,N+1,N+1);
    MultiScalarField3D<T> yCoord(N+1,N+1,N+1);
    MultiScalarField3D<T> zCoord(N+1,N+1,N+1);

    const T pi = (T)4.*std::atan((T)1.);
    setToCoordinate(xCoord, xCoord.getBoundingBox(), 0); multiplyInPlace(xCoord, (T)2.*pi/(T)N);
    setToCoordinate(yCoord, yCoord.getBoundingBox(), 1); multiplyInPlace(yCoord, (T)2.*pi/(T)N);
    setToCoordinate(zCoord, zCoord.getBoundingBox(), 2); multiplyInPlace(zCoord, (T)2.*pi/(T)N);

    MultiScalarField3D<T> wave(N+1,N+1,N+1);

    T(*Tsin)(T) = std::sin;  // Explicit reference to the overloaded "sin" is required for automatic template instantiation.
    T(*Tcos)(T) = std::cos;  // Explicit reference to the overloaded "cos" is required for automatic template instantiation.
    multiply(*evaluate(Tsin, xCoord), *evaluate(Tcos, yCoord), wave, wave.getBoundingBox());

    imageWriter.writeScaledGif("wave", *extractSubDomain(wave, Box3D(0,N, 0,N, N/2,N/2)));
}
