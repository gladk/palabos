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

/** \file
 * A timer class for benchmarking program parts -- header file.
 */
#ifndef PLB_LOGFILES_H
#define PLB_LOGFILES_H

#include <string>
#include <map>

namespace plb {

namespace global {

/// A globally accessible log file.
class PlbLogFile {
public:
    PlbLogFile(std::string fName, bool parallel_);
    ~PlbLogFile();
    /// Start a new section, with corresponding indentation.
    void push(std::string sectionName);
    /// End section, and unindent.
    void pop();
    /// Write a log entry (endline is automatic).
    void entry(std::string entryText);
    /// Write a log entry (endline is automatic) and flush the file buffer.
    void flushEntry(std::string entryText);
private:
    PlbLogFile(PlbLogFile const& rhs);
    PlbLogFile& operator=(PlbLogFile const& rhs);
private:
    bool parallel;
    std::ofstream* ofile;
    int indentation;
    std::string indentSpaces;
private:
friend PlbLogFile& logfile(std::string nameOfLogFile);
friend PlbLogFile& logfile_nonparallel(std::string nameOfLogFile);
};

class LogFileCollection {
public:
    LogFileCollection(bool parallel_);
    ~LogFileCollection();
    PlbLogFile& get(std::string nameOfLogFile);
private:
    LogFileCollection(LogFileCollection const& rhs);
    LogFileCollection& operator=(LogFileCollection const& rhs);
private:
    std::map<std::string, PlbLogFile*> collection;
    bool parallel;
};

PlbLogFile& logfile(std::string nameOfLogFile);
PlbLogFile& logfile_nonparallel(std::string nameOfLogFile);

}  // namespace global

}  // namespace plb

#endif  // PLB_LOGFILES_H
