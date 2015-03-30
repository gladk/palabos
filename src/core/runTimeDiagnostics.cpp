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

#include "core/runTimeDiagnostics.h"
#include "parallelism/mpiManager.h"

namespace plb {

DiagnosticFileSingleton::DiagnosticFileSingleton(std::string fName_)
    : ofile( 0 ),
      fName(fName_)
{ }

DiagnosticFileSingleton::~DiagnosticFileSingleton()
{
    delete ofile;
}

void DiagnosticFileSingleton::write(std::string message) {
    if (global::mpi().isMainProcessor()) {
        if (!ofile) {
            ofile = new std::ofstream( std::string(global::directories().getLogOutDir() + fName).c_str() );
        }
        (*ofile) << message << std::endl;
    }
}

DiagnosticFileSingleton& warningFile() {
    static DiagnosticFileSingleton warningFileSingleton("plbWarning.dat");
    return warningFileSingleton;
}

void plbWarning(bool issueWarning, std::string message) {
#ifdef PLB_MPI_PARALLEL
    int warningCounter = issueWarning ? 1 : 0;
    int warningSum = 0;
    global::mpi().reduce(warningCounter, warningSum, MPI_SUM);
    issueWarning = warningSum > 0;

    if (global::mpi().isMainProcessor()) {
#endif
        if (issueWarning) {
            warningFile().write(message);
        }
#ifdef PLB_MPI_PARALLEL
    }
#endif
}

void plbWarning(std::string message) {
    if (global::mpi().isMainProcessor()) {
        warningFile().write(message);
    }
}


DiagnosticFileSingleton& errorFile() {
    static DiagnosticFileSingleton errorFileSingleton("plbError.dat");
    return errorFileSingleton;
}

void plbMemoryError(bool issueError, std::string message) {
#ifdef PLB_MPI_PARALLEL
    int errorCounter = issueError ? 1 : 0;
    global::mpi().reduceAndBcast(errorCounter,  MPI_SUM);
    issueError = errorCounter > 0;
#endif

    if (issueError) {
        errorFile().write(message);
        throw PlbMemoryException(message);
    }
}

void plbMemoryError(std::string message) {
    throw PlbMemoryException(message);
}

void plbMainProcMemoryError(bool issueError, std::string message) {
    int intIssueError = (int)(issueError);
    global::mpi().bCast(&intIssueError, 1);

    if (intIssueError) {
        throw PlbMemoryException(message);
    }
}


void plbIOError(bool issueError, std::string message) {
#ifdef PLB_MPI_PARALLEL
    int errorCounter = issueError ? 1 : 0;
    global::mpi().reduceAndBcast(errorCounter,  MPI_SUM);
    issueError = errorCounter > 0;
#endif

    if (issueError) {
        errorFile().write(message);
        throw PlbIOException(message);
    }
}

void plbIOError(std::string message) {
    throw PlbIOException(message);
}

void plbMainProcIOError(bool issueError, std::string message) {
    int intIssueError = (int)(issueError);
    global::mpi().bCast(&intIssueError, 1);

    if (intIssueError) {
        throw PlbIOException(message);
    }
}

void plbNetworkError(bool issueError, std::string message) {
#ifdef PLB_MPI_PARALLEL
    int errorCounter = issueError ? 1 : 0;
    global::mpi().reduceAndBcast(errorCounter,  MPI_SUM);
    issueError = errorCounter > 0;
#endif

    if (issueError) {
        errorFile().write(message);
        throw PlbNetworkException(message);
    }
}

void plbNetworkError(std::string message) {
    throw PlbNetworkException(message);
}

void plbMainProcNetworkError(bool issueError, std::string message) {
    int intIssueError = (int)(issueError);
    global::mpi().bCast(&intIssueError, 1);

    if (intIssueError) {
        throw PlbNetworkException(message);
    }
}


void plbLogicError(bool issueError, std::string message) {
#ifdef PLB_MPI_PARALLEL
    int errorCounter = issueError ? 1 : 0;
    global::mpi().reduceAndBcast(errorCounter,  MPI_SUM);
    issueError = errorCounter > 0;
#endif

    if (issueError) {
        errorFile().write(message);
        throw PlbLogicException(message);
    }
}

void plbLogicError(std::string message) {
    throw PlbLogicException(message);
}

void plbMainProcLogicError(bool issueError, std::string message) {
    int intIssueError = (int)(issueError);
    global::mpi().bCast(&intIssueError, 1);

    if (intIssueError) {
        throw PlbLogicException(message);
    }
}


PlbMemoryException::PlbMemoryException(std::string message_) throw()
    : message(std::string("Palabos memory error: ") + message_)
{ }

const char* PlbMemoryException::what() const throw() {

    return message.c_str();
}

PlbIOException::PlbIOException(std::string message_) throw()
    : message(std::string("Palabos IO error: ") + message_)
{ }

const char* PlbIOException::what() const throw() {

    return message.c_str();
}

PlbNetworkException::PlbNetworkException(std::string message_) throw()
    : message(std::string("Palabos network error: ") + message_)
{ }

const char* PlbNetworkException::what() const throw() {

    return message.c_str();
}

PlbLogicException::PlbLogicException(std::string message_) throw()
    : message(std::string("Palabos logic error: ") + message_)
{ }

const char* PlbLogicException::what() const throw() {
    return message.c_str();
}

PlbOutOfRangeException::PlbOutOfRangeException(std::string message_) throw()
    : message(std::string("Palabos out-of-range error: ") + message_)
{ }

const char* PlbOutOfRangeException::what() const throw() {
    return message.c_str();
}

}  // namespace plb
