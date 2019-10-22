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

#include "core/runTimeDiagnostics.h"
#include "core/util.h"
#include "parallelism/mpiManager.h"

#include <cstdio>
#include <csignal>

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

void printSerially(FILE* fp, std::string const& message)
{
    global::mpi().barrier(); // Not strictly needed.
    if (global::mpi().getRank() == 0) {
        fprintf(fp, "%s", message.c_str());
        fflush(fp);
    }
    std::string msg(message);
    for (plint iProc = 1; iProc < global::mpi().getSize(); iProc++) {
        bool iAmRoot = (iProc == global::mpi().getRank());
        global::mpi().sendToMaster(msg, iAmRoot);
        if (global::mpi().getRank() == 0) {
            fprintf(fp, "%s", msg.c_str());
            fflush(fp);
        }
        global::mpi().barrier(); // Not strictly needed.
    }
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


void plbGenericError(bool issueError, std::string message) {
#ifdef PLB_MPI_PARALLEL
    int errorCounter = issueError ? 1 : 0;
    global::mpi().reduceAndBcast(errorCounter, MPI_SUM);
    issueError = errorCounter > 0;
#endif

    if (issueError) {
        errorFile().write(message);
        throw PlbGenericException(message);
    }
}

void plbGenericError(std::string message) {
    throw PlbGenericException(message);
}

void plbMainProcGenericError(bool issueError, std::string message) {
    int intIssueError = (int)(issueError);
    global::mpi().bCast(&intIssueError, 1);

    if (intIssueError) {
        throw PlbGenericException(message);
    }
}


void plbMemoryError(bool issueError, std::string message) {
#ifdef PLB_MPI_PARALLEL
    int errorCounter = issueError ? 1 : 0;
    global::mpi().reduceAndBcast(errorCounter, MPI_SUM);
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
    global::mpi().reduceAndBcast(errorCounter, MPI_SUM);
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
    global::mpi().reduceAndBcast(errorCounter, MPI_SUM);
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
    global::mpi().reduceAndBcast(errorCounter, MPI_SUM);
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


void plbOutOfRangeError(bool issueError, std::string message) {
#ifdef PLB_MPI_PARALLEL
    int errorCounter = issueError ? 1 : 0;
    global::mpi().reduceAndBcast(errorCounter, MPI_SUM);
    issueError = errorCounter > 0;
#endif

    if (issueError) {
        errorFile().write(message);
        throw PlbOutOfRangeException(message);
    }
}

void plbOutOfRangeError(std::string message) {
    throw PlbOutOfRangeException(message);
}

void plbMainProcOutOfRangeError(bool issueError, std::string message) {
    int intIssueError = (int)(issueError);
    global::mpi().bCast(&intIssueError, 1);

    if (intIssueError) {
        throw PlbOutOfRangeException(message);
    }
}


PlbGenericException::PlbGenericException(std::string message_) throw()
    : message(std::string("Palabos generic exception: ") + message_)
{ }

const char* PlbGenericException::what() const throw() {

    return message.c_str();
}


PlbMemoryException::PlbMemoryException(std::string message_) throw()
    : message(std::string("Palabos memory exception: ") + message_)
{ }

const char* PlbMemoryException::what() const throw() {

    return message.c_str();
}

PlbIOException::PlbIOException(std::string message_) throw()
    : message(std::string("Palabos IO exception: ") + message_)
{ }

const char* PlbIOException::what() const throw() {

    return message.c_str();
}

PlbNetworkException::PlbNetworkException(std::string message_) throw()
    : message(std::string("Palabos network exception: ") + message_)
{ }

const char* PlbNetworkException::what() const throw() {

    return message.c_str();
}

PlbLogicException::PlbLogicException(std::string message_) throw()
    : message(std::string("Palabos logic exception: ") + message_)
{ }

const char* PlbLogicException::what() const throw() {
    return message.c_str();
}

PlbOutOfRangeException::PlbOutOfRangeException(std::string message_) throw()
    : message(std::string("Palabos out-of-range exception: ") + message_)
{ }

const char* PlbOutOfRangeException::what() const throw() {
    return message.c_str();
}

namespace global {

void PlbErrors::registerError(std::string const& message)
{
    messages.insert(message);
}

void PlbErrors::registerMemoryError(std::string const& message)
{
    messages.insert(message);
    memoryMessages.insert(message);
}

void PlbErrors::registerIOError(std::string const& message)
{
    messages.insert(message);
    ioMessages.insert(message);
}

void PlbErrors::registerNetworkError(std::string const& message)
{
    messages.insert(message);
    networkMessages.insert(message);
}

void PlbErrors::registerLogicError(std::string const& message)
{
    messages.insert(message);
    logicMessages.insert(message);
}

void PlbErrors::registerOutOfRangeError(std::string const& message)
{
    messages.insert(message);
    outOfRangeMessages.insert(message);
}

void PlbErrors::clear()
{
    messages.clear();
    memoryMessages.clear();
    ioMessages.clear();
    networkMessages.clear();
    logicMessages.clear();
    outOfRangeMessages.clear();
}

void PlbErrors::memoryClear()
{
    for (std::set<std::string>::const_iterator memoryIt = memoryMessages.begin(); memoryIt != memoryMessages.end(); ++memoryIt) {
        std::set<std::string>::iterator it = messages.find(*memoryIt);
        PLB_ASSERT(it != messages.end());
        if (!ioMessageExists(*memoryIt) && !networkMessageExists(*memoryIt) && !logicMessageExists(*memoryIt) &&
                !outOfRangeMessageExists(*memoryIt)) {
            messages.erase(it);
        }
    }
    memoryMessages.clear();
}

void PlbErrors::ioClear()
{
    for (std::set<std::string>::const_iterator ioIt = ioMessages.begin(); ioIt != ioMessages.end(); ++ioIt) {
        std::set<std::string>::iterator it = messages.find(*ioIt);
        PLB_ASSERT(it != messages.end());
        if (!memoryMessageExists(*ioIt) && !networkMessageExists(*ioIt) && !logicMessageExists(*ioIt) &&
                !outOfRangeMessageExists(*ioIt)) {
            messages.erase(it);
        }
    }
    ioMessages.clear();
}

void PlbErrors::networkClear()
{
    for (std::set<std::string>::const_iterator networkIt = networkMessages.begin(); networkIt != networkMessages.end(); ++networkIt) {
        std::set<std::string>::iterator it = messages.find(*networkIt);
        PLB_ASSERT(it != messages.end());
        if (!memoryMessageExists(*networkIt) && !ioMessageExists(*networkIt) && !logicMessageExists(*networkIt) &&
                !outOfRangeMessageExists(*networkIt)) {
            messages.erase(it);
        }
    }
    networkMessages.clear();
}

void PlbErrors::logicClear()
{
    for (std::set<std::string>::const_iterator logicIt = logicMessages.begin(); logicIt != logicMessages.end(); ++logicIt) {
        std::set<std::string>::iterator it = messages.find(*logicIt);
        PLB_ASSERT(it != messages.end());
        if (!memoryMessageExists(*logicIt) && !ioMessageExists(*logicIt) && !networkMessageExists(*logicIt) &&
                !outOfRangeMessageExists(*logicIt)) {
            messages.erase(it);
        }
    }
    logicMessages.clear();
}

void PlbErrors::outOfRangeClear()
{
    for (std::set<std::string>::const_iterator outOfRangeIt = outOfRangeMessages.begin();
            outOfRangeIt != outOfRangeMessages.end(); ++outOfRangeIt) {
        std::set<std::string>::iterator it = messages.find(*outOfRangeIt);
        PLB_ASSERT(it != messages.end());
        if (!memoryMessageExists(*outOfRangeIt) && !ioMessageExists(*outOfRangeIt) && !networkMessageExists(*outOfRangeIt) &&
                !logicMessageExists(*outOfRangeIt)) {
            messages.erase(it);
        }
    }
    outOfRangeMessages.clear();
}

bool PlbErrors::empty() const
{
    return messages.empty();
}

bool PlbErrors::memoryEmpty() const
{
    return memoryMessages.empty();
}

bool PlbErrors::ioEmpty() const
{
    return ioMessages.empty();
}

bool PlbErrors::networkEmpty() const
{
    return networkMessages.empty();
}

bool PlbErrors::logicEmpty() const
{
    return logicMessages.empty();
}

bool PlbErrors::outOfRangeEmpty() const
{
    return outOfRangeMessages.empty();
}

bool PlbErrors::messageExists(std::string const& message) const
{
    return (messages.find(message) != messages.end());
}

bool PlbErrors::memoryMessageExists(std::string const& message) const
{
    return (memoryMessages.find(message) != memoryMessages.end());
}

bool PlbErrors::ioMessageExists(std::string const& message) const
{
    return (ioMessages.find(message) != ioMessages.end());
}

bool PlbErrors::networkMessageExists(std::string const& message) const
{
    return (networkMessages.find(message) != networkMessages.end());
}

bool PlbErrors::logicMessageExists(std::string const& message) const
{
    return (logicMessages.find(message) != logicMessages.end());
}

bool PlbErrors::outOfRangeMessageExists(std::string const& message) const
{
    return (outOfRangeMessages.find(message) != outOfRangeMessages.end());
}

std::string PlbErrors::allMessages() const
{
    std::string message("No error occurred in this process\n");
    if (!empty()) {
        message = "\n  List of error messages (not chronological):";
        std::set<std::string> const& messages = getMessages();
        plint iMessage = 1;
        for (std::set<std::string>::const_iterator it = messages.begin(); it != messages.end(); ++it) {
            message += "\n    " + util::val2str(iMessage) + ") " + *it;
            iMessage++;
        }
        message += "\n";
    }
    return message;
}

std::string PlbErrors::allMemoryMessages() const
{
    std::string message("No memory error occurred in this process\n");
    if (!memoryEmpty()) {
        message = "\n  List of memory error messages (not chronological):";
        std::set<std::string> const& messages = getMemoryMessages();
        plint iMessage = 1;
        for (std::set<std::string>::const_iterator it = messages.begin(); it != messages.end(); ++it) {
            message += "\n    " + util::val2str(iMessage) + ") " + *it;
            iMessage++;
        }
        message += "\n";
    }
    return message;
}

std::string PlbErrors::allIOMessages() const
{
    std::string message("No IO error occurred in this process\n");
    if (!ioEmpty()) {
        message = "\n  List of IO error messages (not chronological):";
        std::set<std::string> const& messages = getIOMessages();
        plint iMessage = 1;
        for (std::set<std::string>::const_iterator it = messages.begin(); it != messages.end(); ++it) {
            message += "\n    " + util::val2str(iMessage) + ") " + *it;
            iMessage++;
        }
        message += "\n";
    }
    return message;
}

std::string PlbErrors::allNetworkMessages() const
{
    std::string message("No network error occurred in this process\n");
    if (!networkEmpty()) {
        message = "\n  List of network error messages (not chronological):";
        std::set<std::string> const& messages = getNetworkMessages();
        plint iMessage = 1;
        for (std::set<std::string>::const_iterator it = messages.begin(); it != messages.end(); ++it) {
            message += "\n    " + util::val2str(iMessage) + ") " + *it;
            iMessage++;
        }
        message += "\n";
    }
    return message;
}

std::string PlbErrors::allLogicMessages() const
{
    std::string message("No logic error occurred in this process\n");
    if (!logicEmpty()) {
        message = "\n  List of logic error messages (not chronological):";
        std::set<std::string> const& messages = getLogicMessages();
        plint iMessage = 1;
        for (std::set<std::string>::const_iterator it = messages.begin(); it != messages.end(); ++it) {
            message += "\n    " + util::val2str(iMessage) + ") " + *it;
            iMessage++;
        }
        message += "\n";
    }
    return message;
}

std::string PlbErrors::allOutOfRangeMessages() const
{
    std::string message("No out-of-range error occurred in this process\n");
    if (!outOfRangeEmpty()) {
        message = "\n  List of out-of-range error messages (not chronological):";
        std::set<std::string> const& messages = getOutOfRangeMessages();
        plint iMessage = 1;
        for (std::set<std::string>::const_iterator it = messages.begin(); it != messages.end(); ++it) {
            message += "\n    " + util::val2str(iMessage) + ") " + *it;
            iMessage++;
        }
        message += "\n";
    }
    return message;
}

void PlbErrors::printMessages() const
{
    if (!empty()) {
        int myRank = (int) global::mpi().getRank();
        fprintf(stderr, "Errors at process %d: %s", myRank, allMessages().c_str()); 
        fflush(stderr); // This is not really needed since stderr is either line buffered or unbuffered.
    }
}

void PlbErrors::printMemoryMessages() const
{
    if (!memoryEmpty()) {
        int myRank = (int) global::mpi().getRank();
        fprintf(stderr, "Memory errors at process %d: %s", myRank, allMemoryMessages().c_str()); 
        fflush(stderr); // This is not really needed since stderr is either line buffered or unbuffered.
    }
}

void PlbErrors::printIOMessages() const
{
    if (!ioEmpty()) {
        int myRank = (int) global::mpi().getRank();
        fprintf(stderr, "IO errors at process %d: %s", myRank, allIOMessages().c_str()); 
        fflush(stderr); // This is not really needed since stderr is either line buffered or unbuffered.
    }
}

void PlbErrors::printNetworkMessages() const
{
    if (!networkEmpty()) {
        int myRank = (int) global::mpi().getRank();
        fprintf(stderr, "Network errors at process %d: %s", myRank, allNetworkMessages().c_str()); 
        fflush(stderr); // This is not really needed since stderr is either line buffered or unbuffered.
    }
}

void PlbErrors::printLogicMessages() const
{
    if (!logicEmpty()) {
        int myRank = (int) global::mpi().getRank();
        fprintf(stderr, "Logic errors at process %d: %s", myRank, allLogicMessages().c_str()); 
        fflush(stderr); // This is not really needed since stderr is either line buffered or unbuffered.
    }
}

void PlbErrors::printOutOfRangeMessages() const
{
    if (!outOfRangeEmpty()) {
        int myRank = (int) global::mpi().getRank();
        fprintf(stderr, "Out-of-range errors at process %d: %s", myRank, allOutOfRangeMessages().c_str()); 
        fflush(stderr); // This is not really needed since stderr is either line buffered or unbuffered.
    }
}

std::set<std::string> const& PlbErrors::getMessages() const
{
    return messages;
}

std::set<std::string> const& PlbErrors::getMemoryMessages() const
{
    return memoryMessages;
}

std::set<std::string> const& PlbErrors::getIOMessages() const
{
    return ioMessages;
}

std::set<std::string> const& PlbErrors::getNetworkMessages() const
{
    return networkMessages;
}

std::set<std::string> const& PlbErrors::getLogicMessages() const
{
    return logicMessages;
}

std::set<std::string> const& PlbErrors::getOutOfRangeMessages() const
{
    return outOfRangeMessages;
}

PlbErrors& plbErrors()
{
    static PlbErrors errors;
    return errors;
}

}  // namespace global

void throwGenericExceptionIfErrorsOccurred()
{
    bool issueError = false;
    int errorCounter = global::plbErrors().empty() ? 0 : 1;
#ifdef PLB_MPI_PARALLEL
    global::mpi().reduceAndBcast(errorCounter, MPI_SUM);
#endif
    issueError = errorCounter > 0;
    if (issueError) {
        throw PlbGenericException(global::plbErrors().allMessages());
    }
}

void throwExceptionIfMemoryErrorsOccurred()
{
    bool issueError = false;
    int errorCounter = global::plbErrors().memoryEmpty() ? 0 : 1;
#ifdef PLB_MPI_PARALLEL
    global::mpi().reduceAndBcast(errorCounter, MPI_SUM);
#endif
    issueError = errorCounter > 0;
    if (issueError) {
        throw PlbMemoryException(global::plbErrors().allMemoryMessages());
    }
}

void throwExceptionIfIOErrorsOccurred()
{
    bool issueError = false;
    int errorCounter = global::plbErrors().ioEmpty() ? 0 : 1;
#ifdef PLB_MPI_PARALLEL
    global::mpi().reduceAndBcast(errorCounter, MPI_SUM);
#endif
    issueError = errorCounter > 0;
    if (issueError) {
        throw PlbIOException(global::plbErrors().allIOMessages());
    }
}

void throwExceptionIfNetworkErrorsOccurred()
{
    bool issueError = false;
    int errorCounter = global::plbErrors().networkEmpty() ? 0 : 1;
#ifdef PLB_MPI_PARALLEL
    global::mpi().reduceAndBcast(errorCounter, MPI_SUM);
#endif
    issueError = errorCounter > 0;
    if (issueError) {
        throw PlbNetworkException(global::plbErrors().allNetworkMessages());
    }
}

void throwExceptionIfLogicErrorsOccurred()
{
    bool issueError = false;
    int errorCounter = global::plbErrors().logicEmpty() ? 0 : 1;
#ifdef PLB_MPI_PARALLEL
    global::mpi().reduceAndBcast(errorCounter, MPI_SUM);
#endif
    issueError = errorCounter > 0;
    if (issueError) {
        throw PlbLogicException(global::plbErrors().allLogicMessages());
    }
}

void throwExceptionIfOutOfRangeErrorsOccurred()
{
    bool issueError = false;
    int errorCounter = global::plbErrors().outOfRangeEmpty() ? 0 : 1;
#ifdef PLB_MPI_PARALLEL
    global::mpi().reduceAndBcast(errorCounter, MPI_SUM);
#endif
    issueError = errorCounter > 0;
    if (issueError) {
        throw PlbOutOfRangeException(global::plbErrors().allOutOfRangeMessages());
    }
}

static void printErrorsSerially(std::string const& message)
{
    global::mpi().barrier();
    if (global::mpi().getRank() == 0) {
        fprintf(stderr, "\nErrors\n"); 
        fprintf(stderr, "------\n"); 
        fflush(stderr); // This is not really needed since stderr is either line buffered or unbuffered.
    }
    printSerially(stderr, message);
}

void abortIfErrorsOccurred()
{
    try {
        throwGenericExceptionIfErrorsOccurred();
    } catch (PlbException& exception) {
        std::string message = "Caught exception in process " + util::val2str(global::mpi().getRank()) + ". " + exception.what() + "\n";
        printErrorsSerially(message);
        exit(-1);
    }
}

void abortIfMemoryErrorsOccurred()
{
    try {
        throwExceptionIfMemoryErrorsOccurred();
    } catch (PlbMemoryException& exception) {
        std::string message = "Caught exception in process " + util::val2str(global::mpi().getRank()) + ". " + exception.what() + "\n";
        printErrorsSerially(message);
        exit(-1);
    }
}

void abortIfIOErrorsOccurred()
{
    try {
        throwExceptionIfIOErrorsOccurred();
    } catch (PlbIOException& exception) {
        std::string message = "Caught exception in process " + util::val2str(global::mpi().getRank()) + ". " + exception.what() + "\n";
        printErrorsSerially(message);
        exit(-1);
    }
}

void abortIfNetworkErrorsOccurred()
{
    try {
        throwExceptionIfNetworkErrorsOccurred();
    } catch (PlbNetworkException& exception) {
        std::string message = "Caught exception in process " + util::val2str(global::mpi().getRank()) + ". " + exception.what() + "\n";
        printErrorsSerially(message);
        exit(-1);
    }
}

void abortIfLogicErrorsOccurred()
{
    try {
        throwExceptionIfLogicErrorsOccurred();
    } catch (PlbLogicException& exception) {
        std::string message = "Caught exception in process " + util::val2str(global::mpi().getRank()) + ". " + exception.what() + "\n";
        printErrorsSerially(message);
        exit(-1);
    }
}

void abortIfOutOfRangeErrorsOccurred()
{
    try {
        throwExceptionIfOutOfRangeErrorsOccurred();
    } catch (PlbOutOfRangeException& exception) {
        std::string message = "Caught exception in process " + util::val2str(global::mpi().getRank()) + ". " + exception.what() + "\n";
        printErrorsSerially(message);
        exit(-1);
    }
}

static void signalHandler(int i)
{
    // Technically in the signal handler we should not manipulate static data and we should use only reentrant functions.
    // In this first implementation of basic signal handling we want to be as simple and portable as possible.
    global::plbErrors().printMessages();
    exit(i);
}

void catchStandardSignals()
{
    (void) signal(SIGABRT, signalHandler);
    (void) signal(SIGFPE, signalHandler);
    (void) signal(SIGILL, signalHandler);
    (void) signal(SIGINT, signalHandler);
    (void) signal(SIGSEGV, signalHandler);
    (void) signal(SIGTERM, signalHandler);
}

}  // namespace plb
