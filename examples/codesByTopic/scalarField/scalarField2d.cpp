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

#include "palabos2D.h"
#include "palabos2D.hh"
#include <cmath>

using namespace plb;
using namespace std;
typedef double T;


int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    const int N = 600;
    ImageWriter<T> imageWriter("leeloo");

    MultiScalarField2D<T> field1(N+1,N+1);
    setToConstant(field1, Box2D(0,N/2, 0,N), (T)1.);
    setToConstant(field1, Box2D(N/2+1,N, 0,N), (T)2.);

    MultiScalarField2D<T> field2(N+1,N+1);
    setToConstant(field2, Box2D(0,N, 0,N/2), (T)1.);
    setToConstant(field2, Box2D(0,N, N/2+1,N), (T)2.);

    MultiScalarField2D<T> field3(N+1,N+1);
    setToCoordinate(field3, field3.getBoundingBox(), 0);

    imageWriter.writeScaledGif("field1", field1);
    imageWriter.writeScaledGif("field2", field2);
    imageWriter.writeScaledGif("field3", field3);
    imageWriter.writeScaledGif("field1+field2", *add(field1, field2));
    imageWriter.writeScaledGif("300*field1+field3", *add(*multiply((T)300.,field1), field3));

    const T pi = (T)4.*std::atan((T)1.);
    T(*Tsin)(T) = std::sin;  // Explicit reference to the overloaded "sin" is required for automatic template instantiation.
    imageWriter.writeScaledGif("sine", *evaluate(Tsin , *multiply((T)2.*pi/(T)N, field3)));
}
