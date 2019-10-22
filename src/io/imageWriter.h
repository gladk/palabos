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

#ifndef IMAGE_WRITER_H
#define IMAGE_WRITER_H

#include "core/globalDefs.h"
#include "atomicBlock/dataField2D.h"
#include "multiBlock/multiDataField2D.h"
#include "atomicBlock/dataField3D.h"
#include "multiBlock/multiDataField3D.h"
#include "io/colormaps.h"
#include <sstream>
#include <iomanip>
#include <vector>

namespace plb {

template<typename T>
class ImageWriter {
public:
    ImageWriter(std::string const& map);
    ImageWriter(std::string const& map, plint colorRange_, plint numColors_);
    void setMap(std::string const& map, plint colorRange_, plint numColors_);

    void writePpm(std::string const& fName,
                  ScalarField2D<T>& field,
                  T minVal, T maxVal) const;
    void writeScaledPpm(std::string const& fName,
                        ScalarField2D<T>& field) const;
    void writeGif(std::string const& fName,
                  ScalarField2D<T>& field,
                  T minVal, T maxVal) const;
    void writeGif(std::string const& fName,
                  ScalarField2D<T>& field,
                  T minVal, T maxVal, plint sizeX, plint sizeY) const;
    void writeScaledGif(std::string const& fName,
                        ScalarField2D<T>& field) const;
    void writeScaledGif(std::string const& fName,
                        ScalarField2D<T>& field,
                        plint sizeX, plint sizeY) const;

    void writePpm(std::string const& fName,
                  MultiScalarField2D<T>& field,
                  T minVal, T maxVal) const;
    void writeScaledPpm(std::string const& fName,
                        MultiScalarField2D<T>& field) const;
    void writeGif(std::string const& fName,
                  MultiScalarField2D<T>& field,
                  T minVal, T maxVal) const;
    void writeGif(std::string const& fName,
                  MultiScalarField2D<T>& field,
                  T minVal, T maxVal, plint sizeX, plint sizeY) const;
    void writeScaledGif(std::string const& fName,
                        MultiScalarField2D<T>& field) const;
    void writeScaledGif(std::string const& fName,
                        MultiScalarField2D<T>& field,
                        plint sizeX, plint sizeY) const;


    void writePpm(std::string const& fName,
                  ScalarField3D<T>& field,
                  T minVal, T maxVal) const;
    void writeScaledPpm(std::string const& fName,
                        ScalarField3D<T>& field) const;
    void writeGif(std::string const& fName,
                  ScalarField3D<T>& field,
                  T minVal, T maxVal) const;
    void writeGif(std::string const& fName,
                  ScalarField3D<T>& field,
                  T minVal, T maxVal, plint sizeX, plint sizeY) const;
    void writeScaledGif(std::string const& fName,
                        ScalarField3D<T>& field) const;
    void writeScaledGif(std::string const& fName,
                        ScalarField3D<T>& field,
                        plint sizeX, plint sizeY) const;

    void writePpm(std::string const& fName,
                  MultiScalarField3D<T>& field,
                  T minVal, T maxVal) const;
    void writeScaledPpm(std::string const& fName,
                        MultiScalarField3D<T>& field) const;
    void writeGif(std::string const& fName,
                  MultiScalarField3D<T>& field,
                  T minVal, T maxVal) const;
    void writeGif(std::string const& fName,
                  MultiScalarField3D<T>& field,
                  T minVal, T maxVal, plint sizeX, plint sizeY) const;
    void writeScaledGif(std::string const& fName,
                        MultiScalarField3D<T>& field) const;
    void writeScaledGif(std::string const& fName,
                        MultiScalarField3D<T>& field,
                        plint sizeX, plint sizeY) const;

private:
    void writePpmImplementation (
        std::string const& fName,
        ScalarField2D<T>& localField, T minVal, T maxVal) const;
    void imageMagickPpmToGif(std::string const& fName) const;
    void imageMagickResize(std::string const& fName,
                           plint sizeX, plint sizeY) const;
private:
    plint colorRange, numColors;
    ColorMap colorMap;
};


////////// Standalone functions ////////////////////////////////////////

inline std::string createFileName(std::string name, plint number, plint width) {
    std::stringstream fNameStream;
    fNameStream << name << std::setfill('0') << std::setw(width) << number;
    return fNameStream.str();
}

}  // namespace plb

#endif
