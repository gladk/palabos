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
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace plb;
using namespace std;


template<typename T>
void plbFileToVtk( std::string fName, std::string identifier, VtkImageOutput2D<T>& vtkOut)
{
    parallelIO::SavedFullMultiBlockSerializer2D *serializer
        = new parallelIO::SavedFullMultiBlockSerializer2D(fName);
    Box2D bbox = serializer->getBoundingBox();
    std::string convertType = serializer->dataType();
    pcout << "Adding the field \"" << identifier << "\" from file " << fName
          << " to the VTK file, with type \"" << convertType << "\""<< std::endl;
    if (convertType=="double") {
        vtkOut.template writeData<double> (
                bbox.getNx(), bbox.getNy(),
                serializer->getCellDim(), serializer, identifier );
    }
    else if (convertType=="float") {
        vtkOut.template writeData<float> (
                bbox.getNx(), bbox.getNy(),
                serializer->getCellDim(), serializer, identifier );
    }
    else if (convertType=="int") {
        vtkOut.template writeData<int> (
                bbox.getNx(), bbox.getNy(),
                serializer->getCellDim(), serializer, identifier );
    }
    else {
        plbIOError("Cannot convert to type "+convertType);
    }
}

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);

    plint numArgs(global::argc());

    std::vector<std::string> fNames;
    try {
        if (numArgs<3) {
            throw(PlbIOException("Too few arguments"));
        }

        for (plint iArg=1; iArg<numArgs; ++iArg) {
            fNames.push_back(global::argv(iArg));
        }
    }
    catch(PlbIOException const& exception) {
        pcout << exception.what() << std::endl;
        pcout << "Syntax: " << (std::string)global::argv(0)
              << " plb_file_name1 [plb_file_name2, plb_file_name3, ...] output_name"
              << std::endl;
        return -1;
    }

    try {
        double dx = 1.;
        VtkImageOutput2D<double> vtkOut(FileName(fNames.back()).getName(), dx);
        for (plint iFile=0; iFile<(plint)fNames.size()-1; ++iFile) {
            plbFileToVtk(fNames[iFile], FileName(fNames[iFile]).getName(), vtkOut);
        }
    }
    catch(PlbIOException const& exception) {
        pcout << exception.what() << std::endl;
    }
}
