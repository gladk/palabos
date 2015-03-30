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

#ifndef RUN_TIME_DIAGNOSTICS_H
#define RUN_TIME_DIAGNOSTICS_H

#include "core/globalDefs.h"
#include <string>
#include <fstream>
#include <exception>

namespace plb {

class DiagnosticFileSingleton {
public:
    DiagnosticFileSingleton(std::string fName_);
    ~DiagnosticFileSingleton();
    void write(std::string message);
private:
    DiagnosticFileSingleton();
    std::ofstream* ofile;
    std::string fName;
};

DiagnosticFileSingleton& warningFile();

/// Invoke this function on all processes. If one process yields true,
///   a warning is issued to the warning file.
void plbWarning(bool issueWarning, std::string message);

/// Invoke this function at least on main process. If the main process
///  yields true, a warning is issued to the warning file.
void plbWarning(std::string message);

/// Invoke this function on all processes. If one process yields true,
///   all processes launch an exception.
void plbMemoryError(bool issueError, std::string message);

/// Invoke this function on all processes. All processes must agree on the
///   fact that there is an error.
void plbMemoryError(std::string message);

/// Invoke this function on all processes. The main process decides if
///   there is an error or not.
void plbMainProcMemoryError(bool issueError, std::string message);

void plbIOError(bool issueError, std::string message);
void plbIOError(std::string message);
void plbMainProcIOError(bool issueError, std::string message);

void plbNetworkError(bool issueError, std::string message);
void plbNetworkError(std::string message);
void plbMainProcNetworkError(bool issueError, std::string message);

void plbLogicError(bool issueError, std::string message);
void plbLogicError(std::string message);
void plbMainProcLogicError(bool issueError, std::string message);

class PlbException : public std::exception
{ };

class PlbMemoryException : public PlbException {
public:
    PlbMemoryException(std::string message_) throw();
    virtual ~PlbMemoryException() throw() { }
    virtual const char* what() const throw();
private:
    std::string message;
};

class PlbIOException : public PlbException {
public:
    PlbIOException(std::string message_) throw();
    virtual ~PlbIOException() throw() { }
    virtual const char* what() const throw();
private:
    std::string message;
};

class PlbNetworkException : public PlbException {
public:
    PlbNetworkException(std::string message_) throw();
    virtual ~PlbNetworkException() throw() { }
    virtual const char* what() const throw();
private:
    std::string message;
};

class PlbLogicException : public PlbException {
public:
    PlbLogicException(std::string message_) throw();
    virtual ~PlbLogicException() throw() { }
    virtual const char* what() const throw();
private:
    std::string message;
};

class PlbOutOfRangeException : public PlbException {
public:
    PlbOutOfRangeException(std::string message_) throw();
    virtual ~PlbOutOfRangeException() throw() { }
    virtual const char* what() const throw();
private:
    std::string message;
};

}  // namespace plb

#endif  // RUN_TIME_DIAGNOSTICS_H
