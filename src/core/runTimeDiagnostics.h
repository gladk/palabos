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

#ifndef RUN_TIME_DIAGNOSTICS_H
#define RUN_TIME_DIAGNOSTICS_H

#include "core/globalDefs.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <exception>
#include <set>

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
DiagnosticFileSingleton& errorFile();

// This function prints messages sent from each process, sorted by increasing
// process rank. All processes must call this function. The messages
// are sent through the master process.
void printSerially(FILE* fp, std::string const& message);

/// Invoke this function on all processes. If one process yields true,
///   a warning is issued to the warning file.
void plbWarning(bool issueWarning, std::string message);

/// Invoke this function at least on main process. If the main process
///  yields true, a warning is issued to the warning file.
void plbWarning(std::string message);

/// Invoke this function on all processes. If one process yields true,
///   all processes launch an exception.
void plbGenericError(bool issueError, std::string message);

/// Invoke this function on all processes. All processes must agree on the
///   fact that there is an error.
void plbGenericError(std::string message);

/// Invoke this function on all processes. The main process decides if
///   there is an error or not.
void plbMainProcGenericError(bool issueError, std::string message);

void plbMemoryError(bool issueError, std::string message);
void plbMemoryError(std::string message);
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

void plbOutOfRangeError(bool issueError, std::string message);
void plbOutOfRangeError(std::string message);
void plbMainProcOutOfRangeError(bool issueError, std::string message);

class PlbException : public std::exception
{ };

class PlbGenericException : public PlbException {
public:
    PlbGenericException(std::string message_) throw();
    virtual ~PlbGenericException() throw() { }
    virtual const char* what() const throw();
private:
    std::string message;
};

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

namespace global {

class PlbErrors {
public:
    void registerError(std::string const& message);
    void registerMemoryError(std::string const& message);
    void registerIOError(std::string const& message);
    void registerNetworkError(std::string const& message);
    void registerLogicError(std::string const& message);
    void registerOutOfRangeError(std::string const& message);

    void clear();
    void memoryClear();
    void ioClear();
    void networkClear();
    void logicClear();
    void outOfRangeClear();

    bool empty() const;
    bool memoryEmpty() const;
    bool ioEmpty() const;
    bool networkEmpty() const;
    bool logicEmpty() const;
    bool outOfRangeEmpty() const;

    bool messageExists(std::string const& message) const;
    bool memoryMessageExists(std::string const& message) const;
    bool ioMessageExists(std::string const& message) const;
    bool networkMessageExists(std::string const& message) const;
    bool logicMessageExists(std::string const& message) const;
    bool outOfRangeMessageExists(std::string const& message) const;

    std::string allMessages() const;
    std::string allMemoryMessages() const;
    std::string allIOMessages() const;
    std::string allNetworkMessages() const;
    std::string allLogicMessages() const;
    std::string allOutOfRangeMessages() const;

    void printMessages() const;
    void printMemoryMessages() const;
    void printIOMessages() const;
    void printNetworkMessages() const;
    void printLogicMessages() const;
    void printOutOfRangeMessages() const;

    std::set<std::string> const& getMessages() const;
    std::set<std::string> const& getMemoryMessages() const;
    std::set<std::string> const& getIOMessages() const;
    std::set<std::string> const& getNetworkMessages() const;
    std::set<std::string> const& getLogicMessages() const;
    std::set<std::string> const& getOutOfRangeMessages() const;
private:
    std::set<std::string> messages;
    std::set<std::string> memoryMessages;
    std::set<std::string> ioMessages;
    std::set<std::string> networkMessages;
    std::set<std::string> logicMessages;
    std::set<std::string> outOfRangeMessages;
};

PlbErrors& plbErrors();

}  // namespace global

// If the exception is caught and handled, and the code does not abort,
// remember to use one of the global::plbErrors().clear() functions afterwards.

void throwGenericExceptionIfErrorsOccurred();
void throwExceptionIfMemoryErrorsOccurred();
void throwExceptionIfIOErrorsOccurred();
void throwExceptionIfNetworkErrorsOccurred();
void throwExceptionIfLogicErrorsOccurred();
void throwExceptionIfOutOfRangeErrorsOccurred();

void abortIfErrorsOccurred();
void abortIfMemoryErrorsOccurred();
void abortIfIOErrorsOccurred();
void abortIfNetworkErrorsOccurred();
void abortIfLogicErrorsOccurred();
void abortIfOutOfRangeErrorsOccurred();

void catchStandardSignals();

}  // namespace plb

#endif  // RUN_TIME_DIAGNOSTICS_H
