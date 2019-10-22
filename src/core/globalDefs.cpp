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

#include "core/globalDefs.h"

namespace plb {

namespace global {

IOpolicyClass::IOpolicyClass()
    : streamOrdering(IndexOrdering::forward),
      endianSwitchOnBase64out(false),
      endianSwitchOnBase64in(false),
      stlLowerBoundFlag(false),
      stlLowerBound(-1.),
      parallelIOflag(true)
{ }

void IOpolicyClass::setIndexOrderingForStreams(IndexOrdering::OrderingT streamOrdering_) {
    streamOrdering = streamOrdering_;
}

IndexOrdering::OrderingT IOpolicyClass::getIndexOrderingForStreams() const {
    return streamOrdering;
}

void IOpolicyClass::setEndianSwitchOnBase64out(bool doSwitch) {
    endianSwitchOnBase64out = doSwitch;
}

bool IOpolicyClass::getEndianSwitchOnBase64out() {
    return endianSwitchOnBase64out;
}

void IOpolicyClass::setEndianSwitchOnBase64in(bool doSwitch) {
    endianSwitchOnBase64in = doSwitch;
}

bool IOpolicyClass::getEndianSwitchOnBase64in() {
    return endianSwitchOnBase64in;
}

void IOpolicyClass::setStlFilesHaveLowerBound(bool flag) {
    stlLowerBoundFlag = flag;
}

bool IOpolicyClass::stlFilesHaveLowerBound() const {
    return stlLowerBoundFlag;
}

void IOpolicyClass::setLowerBoundForStlFiles(double stlLowerBound_) {
    stlLowerBound = stlLowerBound_;
}

double IOpolicyClass::getLowerBoundForStlFiles() const {
    return stlLowerBound;
}

void IOpolicyClass::activateParallelIO(bool activate) {
    parallelIOflag = activate;
}

bool IOpolicyClass::useParallelIO() const {
    return parallelIOflag;
}

/** Directories are default initialized to working directory.
 */
Directories::Directories()
{
    setOutputDir("");
    setInputDir("");
}

void Directories::setOutputDir(std::string outputDir_) {
    outputDir = outputDir_;
    setLogOutDir(outputDir);
    setImageOutDir(outputDir);
    setVtkOutDir(outputDir);
}

void Directories::setLogOutDir(std::string logOutDir_) {
    logOutDir = logOutDir_;
}

void Directories::setImageOutDir(std::string imageOutDir_) {
    imageOutDir = imageOutDir_;
}

void Directories::setVtkOutDir(std::string vtkOutDir_) {
    vtkOutDir = vtkOutDir_;
}

/** The value of input dir is not used automatically inside
 *  Palabos, but you can refer to it explicitly when reading
 *  input files.
 */
void Directories::setInputDir(std::string inputDir_) {
    inputDir = inputDir_;
}

std::string Directories::getLogOutDir() const {
    return logOutDir;
}

std::string Directories::getImageOutDir() const {
    return imageOutDir;
}

std::string Directories::getVtkOutDir() const {
    return vtkOutDir;
}

std::string Directories::getInputDir() const {
    return inputDir;
}

std::string Directories::getOutputDir() const {
    return outputDir;
}

}  // namespace global

}  // namespace plb
