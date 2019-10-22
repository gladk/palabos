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

#include "parallelism/mpiManager.h"
#include "core/plbLogFiles.h"
#include "core/util.h"
#include <utility>
#include <fstream>

namespace plb {

namespace global {

PlbLogFile::PlbLogFile(std::string fName, bool parallel_)
    : parallel(parallel_),
      ofile(0),
      indentation(0),
      indentSpaces("")
{
    if (parallel) {
        if (global::mpi().isMainProcessor()) {
            ofile = new std::ofstream (
                    (global::directories().getLogOutDir()+fName ).c_str() );
        }
    }
    else {
        ofile = new std::ofstream (
                ( global::directories().getLogOutDir() +
                  util::val2str(global::mpi().getRank())+"_"+fName ).c_str() );
    }
}

PlbLogFile::~PlbLogFile() {
    delete ofile;
}

PlbLogFile::PlbLogFile(PlbLogFile const& rhs)
{ }

PlbLogFile& PlbLogFile::operator=(PlbLogFile const& rhs)
{
    return *this;
}

void PlbLogFile::push(std::string sectionName)
{
    entry("\nSECTION " + sectionName);
    indentation += 4;
    indentSpaces = std::string(indentation, ' ');
}

void PlbLogFile::pop()
{
    indentation -= 4;
    indentSpaces = std::string(indentation, ' ');
}

void PlbLogFile::entry(std::string entryText) {
    if (ofile) {
        (*ofile) << indentSpaces << entryText << "\n";
    }
}

void PlbLogFile::flushEntry(std::string entryText) {
    if (ofile) {
        // Using std::endl enforces file-buffer flush.
        (*ofile) << indentSpaces << entryText << std::endl;
    }
}

LogFileCollection::LogFileCollection(bool parallel_)
    : parallel(parallel_)
{ }

LogFileCollection::~LogFileCollection() {
    std::map<std::string, PlbLogFile*>::iterator it=collection.begin();
    for (; it != collection.end(); ++it) {
        delete it->second;
    }
}

LogFileCollection::LogFileCollection(LogFileCollection const& rhs)
{ }

LogFileCollection& LogFileCollection::operator=(LogFileCollection const& rhs) {
    return *this;
}

PlbLogFile& LogFileCollection::get(std::string nameOfLogFile) {
    std::map<std::string, PlbLogFile*>::iterator it=collection.find(nameOfLogFile);
    if (it == collection.end()) {
        PlbLogFile* logfile = new PlbLogFile(nameOfLogFile, parallel);
        collection.insert(make_pair(nameOfLogFile, logfile));
        return *logfile;
    }
    else {
        return *it->second;
    }
}

PlbLogFile& logfile(std::string nameOfLogFile) {
    static LogFileCollection collection(true);
    return collection.get(nameOfLogFile);
}

PlbLogFile& logfile_nonparallel(std::string nameOfLogFile) {
    static LogFileCollection collection(false);
    return collection.get(nameOfLogFile);
}

}  // namespace global

}  // namespace plb
