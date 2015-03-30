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

/** \file
 * Wrapper functions that simplify the use of MPI, template instatiations
 */

#ifdef PLB_MPI_PARALLEL

#include "parallelism/mpiManager.h"
#include "core/plbDebug.h"
#include "core/plbComplex.h"
#include "core/plbComplex.hh"
#include <algorithm>
#include <iostream>

namespace plb {

namespace global {

MpiManager::MpiManager()
    : ok(false),
      responsibleForMpiMachine(false)
{ }

MpiManager::~MpiManager() {
    if (responsibleForMpiMachine) {
        MPI_Finalize();
        ok = false;
        responsibleForMpiMachine = false;
    }
}

void MpiManager::init(int *argc, char ***argv, bool verbous) {
    if (verbous) {
        std::cerr << "Constructing an MPI thread" << std::endl;
    }
    int ok1 = MPI_Init(argc, argv);
    // If I'm the one who calls MPI_Init, then I need to be
    // the one who calls MPI_Finalize.
    responsibleForMpiMachine = true;
    globalCommunicator = MPI_COMM_WORLD;
    int ok2 = MPI_Comm_rank(getGlobalCommunicator(),&taskId);
    int ok3 = MPI_Comm_size(getGlobalCommunicator(),&numTasks);
    ok = (ok1==0 && ok2==0 && ok3==0);
}

void MpiManager::init(MPI_Comm globalCommunicator_) {
    globalCommunicator = globalCommunicator_;
    int ok1 = MPI_Comm_rank(getGlobalCommunicator(),&taskId);
    int ok2 = MPI_Comm_size(getGlobalCommunicator(),&numTasks);
    ok = (ok1==0 && ok2==0);
}

void MpiManager::init() {
    init(MPI_COMM_WORLD);
}

int MpiManager::getSize() const {
    return numTasks;
}

int MpiManager::getRank() const {
    return taskId;
}

int MpiManager::bossId() const {
    return 0;
}

bool MpiManager::isMainProcessor() const {
    return bossId() == getRank();
}

double MpiManager::getTime() const {
    if (!ok) return 0.;
    return MPI_Wtime();
}

MPI_Comm MpiManager::getGlobalCommunicator() const {
    return globalCommunicator;
}

void MpiManager::barrier() {
    if (!ok) return;
    MPI_Barrier(getGlobalCommunicator());
}

template <>
void MpiManager::send<char>(char *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::send<int>(int *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_INT, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::send<bool>(bool *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count*sizeof(bool), MPI_CHAR, dest, tag, getGlobalCommunicator());
}

#ifdef PLB_BGP
template <>
void MpiManager::send<long long>(long long *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count*sizeof(long long), MPI_CHAR, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::send<unsigned long long>(unsigned long long *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count*sizeof(unsigned long long), MPI_CHAR, dest, tag, getGlobalCommunicator());
}
#endif

template <>
void MpiManager::send<long>(long *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_LONG, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::send<float>(float *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::send<double>(double *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::send<long double>(long double *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_LONG_DOUBLE, dest, tag, getGlobalCommunicator());
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::send<__float128 >(__float128  *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count*sizeof(__float128), MPI_CHAR, dest, tag, getGlobalCommunicator());
}
#endif

template <>
void MpiManager::send<Complex<double> >(Complex<double>  *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count*sizeof(Complex<double>), MPI_CHAR, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::send<Complex<long double> >(Complex<long double>  *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count*sizeof(Complex<long double>), MPI_CHAR, dest, tag, getGlobalCommunicator());
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::send<Complex<__float128> >(Complex<__float128>  *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count*sizeof(Complex<__float128>), MPI_CHAR, dest, tag, getGlobalCommunicator());
}
#endif

template <>
void MpiManager::iSend<char>
    (char *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iSend<int>
    (int *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_INT, dest, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iSend<bool>
    (bool *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(bool), MPI_CHAR, dest, tag, getGlobalCommunicator(), request);
    }
}

#ifdef PLB_BGP
template <>
void MpiManager::iSend<long long>
    (long long *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(long long), MPI_CHAR, dest, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iSend<unsigned long long>
    (unsigned long long *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(unsigned long long), MPI_CHAR, dest, tag, getGlobalCommunicator(), request);
    }
}
#endif

template <>
void MpiManager::iSend<long>
    (long *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_LONG, dest, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iSend<float>
    (float *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iSend<double>
    (double *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iSend<long double>
    (long double *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_LONG_DOUBLE, dest, tag, getGlobalCommunicator(), request);
    }
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::iSend<__float128 >
    (__float128  *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(__float128), MPI_CHAR, dest, tag, getGlobalCommunicator(), request);
    }
}
#endif

template <>
void MpiManager::iSend<Complex<double> >
    (Complex<double>  *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(Complex<double>), MPI_CHAR, dest, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iSend<Complex<long double> >
    (Complex<long double>  *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(Complex<long double>), MPI_CHAR, dest, tag, getGlobalCommunicator(), request);
    }
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::iSend<Complex<__float128> >
    (Complex<__float128>  *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(Complex<__float128>), MPI_CHAR, dest, tag, getGlobalCommunicator(), request);
    }
}
#endif

template <>
void MpiManager::rSend<char>(char *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::rSend<int>(int *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count, MPI_INT, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::rSend<bool>(bool *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count*sizeof(bool), MPI_CHAR, dest, tag, getGlobalCommunicator());
}

#ifdef PLB_BGP
template <>
void MpiManager::rSend<long long>(long long *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count*sizeof(long long), MPI_CHAR, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::rSend<unsigned long long>(unsigned long long *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count*sizeof(unsigned long long), MPI_CHAR, dest, tag, getGlobalCommunicator());
}
#endif

template <>
void MpiManager::rSend<long>(long *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count, MPI_LONG, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::rSend<float>(float *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::rSend<double>(double *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::rSend<long double>(long double *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count, MPI_LONG_DOUBLE, dest, tag, getGlobalCommunicator());
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::rSend<__float128 >(__float128  *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count*sizeof(__float128), MPI_CHAR, dest, tag, getGlobalCommunicator());
}
#endif

template <>
void MpiManager::rSend<Complex<double> >(Complex<double>  *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count*sizeof(Complex<double>), MPI_CHAR, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::rSend<Complex<long double> >(Complex<long double>  *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count*sizeof(Complex<long double>), MPI_CHAR, dest, tag, getGlobalCommunicator());
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::rSend<Complex<__float128> >(Complex<__float128>  *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count*sizeof(Complex<__float128>), MPI_CHAR, dest, tag, getGlobalCommunicator());
}
#endif

template <>
void MpiManager::iSendRequestFree<char>
    (char *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<int>
    (int *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_INT, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<bool>
    (bool *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(bool), MPI_CHAR, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

#ifdef PLB_BGP
template <>
void MpiManager::iSendRequestFree<long long>
    (long long *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(long long), MPI_CHAR, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<unsigned long long>
    (unsigned long long *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(unsigned long long), MPI_CHAR, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}
#endif

template <>
void MpiManager::iSendRequestFree<long>
    (long *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_LONG, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<float>
    (float *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<double>
    (double *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<long double>
    (long double *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_LONG_DOUBLE, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::iSendRequestFree<__float128 >
    (__float128  *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(__float128), MPI_CHAR, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}
#endif

template <>
void MpiManager::iSendRequestFree<Complex<double> >
    (Complex<double>  *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(Complex<double>), MPI_CHAR, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<Complex<long double> >
    (Complex<long double>  *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(Complex<long double>), MPI_CHAR, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::iSendRequestFree<Complex<__float128> >
    (Complex<__float128>  *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(Complex<__float128>), MPI_CHAR, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}
#endif

template <>
void MpiManager::receive<char>(char *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_CHAR, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::receive<int>(int *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_INT, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::receive<bool>(bool *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count*sizeof(bool), MPI_CHAR, source, tag, getGlobalCommunicator(), &status);
}

#ifdef PLB_BGP
template <>
void MpiManager::receive<long long>(long long *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count*sizeof(long long), MPI_CHAR, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::receive<unsigned long long>(unsigned long long *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count*sizeof(unsigned long long), MPI_CHAR, source, tag, getGlobalCommunicator(), &status);
}
#endif

template <>
void MpiManager::receive<long>(long *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_LONG, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::receive<float>(float *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_FLOAT, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::receive<double>(double *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_DOUBLE, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::receive<long double>(long double *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_LONG_DOUBLE, source, tag, getGlobalCommunicator(), &status);
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::receive<__float128 >(__float128 *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count*sizeof(__float128), MPI_CHAR, source, tag, getGlobalCommunicator(), &status);
}
#endif

template <>
void MpiManager::receive<Complex<double> >(Complex<double> *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count*sizeof(Complex<double>), MPI_CHAR, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::receive<Complex<long double> >(Complex<long double> *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count*sizeof(Complex<long double>), MPI_CHAR, source, tag, getGlobalCommunicator(), &status);
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::receive<Complex<__float128> >(Complex<__float128> *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count*sizeof(Complex<__float128>), MPI_CHAR, source, tag, getGlobalCommunicator(), &status);
}
#endif

template <>
void MpiManager::sendToMaster<char>(char* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<int>(int* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<bool>(bool* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

#ifdef PLB_BGP
template <>
void MpiManager::sendToMaster<long long>(long long* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<unsigned long long>(unsigned long long* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}
#endif

template <>
void MpiManager::sendToMaster<long>(long* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<float>(float* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<double>(double* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<long double>(long double* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::sendToMaster<__float128 >(__float128* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}
#endif

template <>
void MpiManager::sendToMaster<Complex<double> >(Complex<double>* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<Complex<long double> >(Complex<long double>* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::sendToMaster<Complex<__float128> >(Complex<__float128>* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}
#endif


template <>
void MpiManager::iRecv<char>(char *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_CHAR, source, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iRecv<int>(int *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_INT, source, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iRecv<bool>(bool *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count*sizeof(bool), MPI_CHAR, source, tag, getGlobalCommunicator(), request);
    }
}

#ifdef PLB_BGP
template <>
void MpiManager::iRecv<long long>(long long *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count*sizeof(long long), MPI_CHAR, source, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iRecv<unsigned long long>(unsigned long long *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count*sizeof(unsigned long long), MPI_CHAR, source, tag, getGlobalCommunicator(), request);
    }
}
#endif

template <>
void MpiManager::iRecv<long>(long *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_LONG, source, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iRecv<float>(float *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_FLOAT, source, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iRecv<double>(double *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_DOUBLE, source, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iRecv<long double>(long double *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_DOUBLE, source, tag, getGlobalCommunicator(), request);
    }
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::iRecv<__float128 >(__float128 *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count*sizeof(__float128), MPI_CHAR, source, tag, getGlobalCommunicator(), request);
    }
}
#endif

template <>
void MpiManager::iRecv<Complex<double> >(Complex<double> *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count*sizeof(Complex<double>), MPI_CHAR, source, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iRecv<Complex<long double> >(Complex<long double> *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count*sizeof(Complex<long double>), MPI_CHAR, source, tag, getGlobalCommunicator(), request);
    }
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::iRecv<Complex<__float128> >(Complex<__float128> *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count*sizeof(Complex<__float128>), MPI_CHAR, source, tag, getGlobalCommunicator(), request);
    }
}
#endif

template <>
void MpiManager::sendRecv<char>
    (char *sendBuf, char *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_CHAR, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_CHAR, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::sendRecv<int>
    (int *sendBuf, int *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_INT, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_INT, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::sendRecv<bool>
    (bool *sendBuf, bool *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count*sizeof(bool),
                 MPI_CHAR, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_LONG, source, tag, getGlobalCommunicator(), &status);
}

#ifdef PLB_BGP
template <>
void MpiManager::sendRecv<long long>
    (long long *sendBuf, long long *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count*sizeof(long long),
                 MPI_CHAR, dest, tag,
                 static_cast<void*>(recvBuf),
                 count*sizeof(long long),
                 MPI_CHAR, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::sendRecv<unsigned long long>
    (unsigned long long *sendBuf, unsigned long long *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count*sizeof(unsigned long long),
                 MPI_CHAR, dest, tag,
                 static_cast<void*>(recvBuf),
                 count*sizeof(long long),
                 MPI_CHAR, source, tag, getGlobalCommunicator(), &status);
}
#endif

template <>
void MpiManager::sendRecv<long>
    (long *sendBuf, long *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_LONG, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_LONG, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::sendRecv<float>
    (float *sendBuf, float *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_FLOAT, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_FLOAT, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::sendRecv<double>
    (double *sendBuf, double *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_DOUBLE, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_DOUBLE, source, tag, getGlobalCommunicator(), &status);

}

template <>
void MpiManager::sendRecv<long double>
    (long double *sendBuf, long double *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_LONG_DOUBLE, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_LONG_DOUBLE, source, tag, getGlobalCommunicator(), &status);

}

#ifdef PLB_USE_GCC
template <>
void MpiManager::sendRecv<__float128 >
    (__float128 *sendBuf, __float128 *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count*sizeof(__float128),
                 MPI_CHAR, dest, tag,
                 static_cast<void*>(recvBuf),
                 count*sizeof(__float128),
                 MPI_CHAR, source, tag, getGlobalCommunicator(), &status);

}
#endif

template <>
void MpiManager::sendRecv<Complex<double> >
    (Complex<double> *sendBuf, Complex<double> *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count*sizeof(Complex<double>),
                 MPI_CHAR, dest, tag,
                 static_cast<void*>(recvBuf),
                 count*sizeof(Complex<double>),
                 MPI_CHAR, source, tag, getGlobalCommunicator(), &status);

}

template <>
void MpiManager::sendRecv<Complex<long double> >
    (Complex<long double> *sendBuf, Complex<long double> *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count*sizeof(Complex<long double>),
                 MPI_CHAR, dest, tag,
                 static_cast<void*>(recvBuf),
                 count*sizeof(Complex<long double>),
                 MPI_CHAR, source, tag, getGlobalCommunicator(), &status);

}

#ifdef PLB_USE_GCC
template <>
void MpiManager::sendRecv<Complex<__float128> >
    (Complex<__float128> *sendBuf, Complex<__float128> *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count*sizeof(Complex<__float128>),
                 MPI_CHAR, dest, tag,
                 static_cast<void*>(recvBuf),
                 count*sizeof(Complex<__float128>),
                 MPI_CHAR, source, tag, getGlobalCommunicator(), &status);

}
#endif

template <>
void MpiManager::scatterv_impl<char>(char* sendBuf, int* sendCounts, int* displs,
                                     char* recvBuf, int recvCount, int root)
{
    if (!ok) return;
    MPI_Scatterv(static_cast<void*>(sendBuf),
                 sendCounts, displs, MPI_CHAR,
                 static_cast<void*>(recvBuf),
                 recvCount, MPI_CHAR, root, getGlobalCommunicator());
}

template <>
void MpiManager::scatterv_impl<int>(int *sendBuf, int* sendCounts, int* displs,
                                int* recvBuf, int recvCount, int root)
{
    if (!ok) return;
    MPI_Scatterv(static_cast<void*>(sendBuf),
                 sendCounts, displs, MPI_INT,
                 static_cast<void*>(recvBuf),
                 recvCount, MPI_INT, root, getGlobalCommunicator());
}

template <>
void MpiManager::scatterv_impl<long>(long *sendBuf, int* sendCounts, int* displs,
                                     long* recvBuf, int recvCount, int root)
{
    if (!ok) return;
    MPI_Scatterv(static_cast<void*>(sendBuf),
                 sendCounts, displs, MPI_LONG,
                 static_cast<void*>(recvBuf),
                 recvCount, MPI_LONG, root, getGlobalCommunicator());
}

template <>
void MpiManager::scatterv_impl<float>(float *sendBuf, int* sendCounts, int* displs,
                                      float* recvBuf, int recvCount, int root)
{
    if (!ok) return;
    MPI_Scatterv(static_cast<void*>(sendBuf),
                 sendCounts, displs, MPI_FLOAT,
                 static_cast<void*>(recvBuf),
                 recvCount, MPI_FLOAT, root, getGlobalCommunicator());
}

template <>
void MpiManager::scatterv_impl<double>(double *sendBuf, int* sendCounts, int* displs,
                                       double* recvBuf, int recvCount, int root)
{
    if (!ok) return;
    MPI_Scatterv(static_cast<void*>(sendBuf),
                 sendCounts, displs, MPI_DOUBLE,
                 static_cast<void*>(recvBuf),
                 recvCount, MPI_DOUBLE, root, getGlobalCommunicator());
}

template <>
void MpiManager::scatterv_impl<long double>(long double *sendBuf, int* sendCounts, int* displs,
                                            long double* recvBuf, int recvCount, int root)
{
    if (!ok) return;
    MPI_Scatterv(static_cast<void*>(sendBuf),
                 sendCounts, displs, MPI_LONG_DOUBLE,
                 static_cast<void*>(recvBuf),
                 recvCount, MPI_LONG_DOUBLE, root, getGlobalCommunicator());
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::scatterv_impl<__float128 >(__float128 *sendBuf, int* sendCounts, int* displs,
                                       __float128* recvBuf, int recvCount, int root)
{
    if (!ok) return;
    PLB_ASSERT( false ); // Not yet implemented.
}
#endif

template <>
void MpiManager::scatterv_impl<Complex<double> >(Complex<double> *sendBuf, int* sendCounts, int* displs,
                                       Complex<double>* recvBuf, int recvCount, int root)
{
    if (!ok) return;
    PLB_ASSERT( false ); // Not yet implemented.
}

template <>
void MpiManager::scatterv_impl<Complex<long double> >(Complex<long double> *sendBuf, int* sendCounts, int* displs,
                                                      Complex<long double>* recvBuf, int recvCount, int root)
{
    if (!ok) return;
    PLB_ASSERT( false ); // Not yet implemented.
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::scatterv_impl<Complex<__float128> >(Complex<__float128> *sendBuf, int* sendCounts, int* displs,
                                       Complex<__float128>* recvBuf, int recvCount, int root)
{
    if (!ok) return;
    PLB_ASSERT( false ); // Not yet implemented.
}
#endif

template <>
void MpiManager::gatherv_impl<char>(char* sendBuf, int sendCount,
                                    char* recvBuf, int* recvCounts, int* displs,
                                    int root)
{
    if (!ok) return;
    MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_CHAR,
                static_cast<void*>(recvBuf), recvCounts, displs, MPI_CHAR,
                root, getGlobalCommunicator());
}

template <>
void MpiManager::gatherv_impl<int>(int* sendBuf, int sendCount,
                                   int* recvBuf, int* recvCounts, int* displs,
                                   int root)
{
    if (!ok) return;
    MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_INT,
                static_cast<void*>(recvBuf), recvCounts, displs, MPI_INT,
                root, getGlobalCommunicator());
}

template <>
void MpiManager::gatherv_impl<long>(long* sendBuf, int sendCount,
                                    long* recvBuf, int* recvCounts, int* displs,
                                    int root)
{
    if (!ok) return;
    MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_LONG,
                static_cast<void*>(recvBuf), recvCounts, displs, MPI_LONG,
                root, getGlobalCommunicator());
}

template <>
void MpiManager::gatherv_impl<float>(float* sendBuf, int sendCount,
                                 float* recvBuf, int* recvCounts, int* displs,
                                 int root)
{
    if (!ok) return;
    MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_FLOAT,
                static_cast<void*>(recvBuf), recvCounts, displs, MPI_FLOAT,
                root, getGlobalCommunicator());
}

template <>
void MpiManager::gatherv_impl<double>(double* sendBuf, int sendCount,
                                      double* recvBuf, int* recvCounts, int* displs,
                                      int root)
{
    if (!ok) return;
    MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_DOUBLE,
                static_cast<void*>(recvBuf), recvCounts, displs, MPI_DOUBLE,
                root, getGlobalCommunicator());
}

template <>
void MpiManager::gatherv_impl<long double>(long double* sendBuf, int sendCount,
                                           long double* recvBuf, int* recvCounts, int* displs,
                                           int root)
{
    if (!ok) return;
    MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_LONG_DOUBLE,
                static_cast<void*>(recvBuf), recvCounts, displs, MPI_LONG_DOUBLE,
                root, getGlobalCommunicator());
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::gatherv_impl<__float128 >(__float128 * sendBuf, int sendCount,
                                  __float128 * recvBuf, int* recvCounts, int* displs,
                                  int root)
{
    if (!ok) return;
    PLB_ASSERT( false ); // Not yet implemented.
}
#endif

template <>
void MpiManager::gatherv_impl<Complex<double> >(Complex<double> * sendBuf, int sendCount,
                                  Complex<double> * recvBuf, int* recvCounts, int* displs,
                                  int root)
{
    if (!ok) return;
    PLB_ASSERT( false ); // Not yet implemented.
}

template <>
void MpiManager::gatherv_impl<Complex<long double> >(Complex<long double> * sendBuf, int sendCount,
                                  Complex<long double> * recvBuf, int* recvCounts, int* displs,
                                  int root)
{
    if (!ok) return;
    PLB_ASSERT( false ); // Not yet implemented.
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::gatherv_impl<Complex<__float128> >(Complex<__float128> * sendBuf, int sendCount,
                                  Complex<__float128> * recvBuf, int* recvCounts, int* displs,
                                  int root)
{
    if (!ok) return;
    PLB_ASSERT( false ); // Not yet implemented.
}
#endif

template <>
void MpiManager::bCast<char>(char* sendBuf, int sendCount, int root)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount, MPI_CHAR, root, getGlobalCommunicator());
}

template <> void MpiManager::bCast<int>(int* sendBuf, int sendCount, int root)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount, MPI_INT, root, getGlobalCommunicator());
}

void MpiManager::bCast(std::string& message, int root)
{
    if (!ok) return;
    int length = (int) message.size();
    bCast(&length, 1, root);
    char* buffer = new char[length+1];
    if (getRank()==root) {
        std::copy(message.c_str(), message.c_str()+length+1, buffer);
    }
    bCast(buffer, length+1, root);
    if (getRank()!=root) {
        message = buffer;
    }
    delete [] buffer;
}

template <>
void MpiManager::bCast<bool>(bool* sendBuf, int sendCount, int root)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount*sizeof(bool), MPI_CHAR, root, getGlobalCommunicator());
}

#ifdef PLB_BGP
template <>
void MpiManager::bCast<long long>(long long* sendBuf, int sendCount, int root)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount*sizeof(long long), MPI_CHAR, root, getGlobalCommunicator());
}

template <>
void MpiManager::bCast<unsigned long long>(unsigned long long* sendBuf, int sendCount, int root)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount*sizeof(unsigned long long), MPI_CHAR, root, getGlobalCommunicator());
}
#endif

template <>
void MpiManager::bCast<long>(long* sendBuf, int sendCount, int root)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount, MPI_LONG, root, getGlobalCommunicator());
}

template <>
void MpiManager::bCast<float>(float* sendBuf, int sendCount, int root)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount, MPI_FLOAT, root, getGlobalCommunicator());
}

template <>
void MpiManager::bCast<double>(double* sendBuf, int sendCount, int root)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount, MPI_DOUBLE, root, getGlobalCommunicator());
}

template <>
void MpiManager::bCast<long double>(long double* sendBuf, int sendCount, int root)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount, MPI_LONG_DOUBLE, root, getGlobalCommunicator());
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::bCast<__float128 >(__float128* sendBuf, int sendCount, int root)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount*sizeof(__float128), MPI_CHAR, root, getGlobalCommunicator());
}
#endif

template <>
void MpiManager::bCast<Complex<double> >(Complex<double>* sendBuf, int sendCount, int root)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount*sizeof(Complex<double>), MPI_CHAR, root, getGlobalCommunicator());
}

template <>
void MpiManager::bCast<Complex<long double> >(Complex<long double>* sendBuf, int sendCount, int root)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount*sizeof(Complex<long double>), MPI_CHAR, root, getGlobalCommunicator());
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::bCast<Complex<__float128> >(Complex<__float128>* sendBuf, int sendCount, int root)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount*sizeof(Complex<__float128>), MPI_CHAR, root, getGlobalCommunicator());
}
#endif

template <>
void MpiManager::bCastThroughMaster<char>(char* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<int>(int* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<bool>(bool* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

#ifdef PLB_BGP
template <>
void MpiManager::bCastThroughMaster<long long>(long long* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<unsigned long long>(unsigned long long* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}
#endif

template <>
void MpiManager::bCastThroughMaster<long>(long* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<float>(float* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<double>(double* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<long double>(long double* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::bCastThroughMaster<__float128 >(__float128 * sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}
#endif

template <>
void MpiManager::bCastThroughMaster<Complex<double> >(Complex<double> * sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<Complex<long double> >(Complex<long double> * sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::bCastThroughMaster<Complex<__float128> >(Complex<__float128> * sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}
#endif

template <>
void MpiManager::reduce<char>(char sendVal, char& recvVal,  MPI_Op op, int root)
{
    if (!ok) return;
    int tmpSend = sendVal;
    int tmpRecv;
    MPI_Reduce(static_cast<void*>(&tmpSend),
               static_cast<void*>(&tmpRecv), 1, MPI_INT, op, root, getGlobalCommunicator());
    recvVal = tmpRecv;
}

template <>
void MpiManager::reduce<int>(int sendVal, int& recvVal,  MPI_Op op, int root)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 1, MPI_INT, op, root, getGlobalCommunicator());
}

template <>
void MpiManager::reduce<long>(long sendVal, long& recvVal,  MPI_Op op, int root)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 1, MPI_LONG, op, root, getGlobalCommunicator());
}

#ifdef PLB_BGP
template <>
void MpiManager::reduce<long long>(long long sendVal, long long& recvVal,  MPI_Op op, int root)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 1, MPI_LONG_LONG, op, root, getGlobalCommunicator());
}

template <>
void MpiManager::reduce<unsigned long long>(unsigned long long sendVal, unsigned long long& recvVal,  MPI_Op op, int root)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 1, MPI_UNSIGNED_LONG_LONG, op, root, getGlobalCommunicator());
}
#endif

template <>
void MpiManager::reduce<float>(float sendVal, float& recvVal,  MPI_Op op, int root)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 1, MPI_FLOAT, op, root, getGlobalCommunicator());
}

template <>
void MpiManager::reduce<double>(double sendVal, double& recvVal,  MPI_Op op, int root)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 1, MPI_DOUBLE, op, root, getGlobalCommunicator());
}

template <>
void MpiManager::reduce<long double>(long double sendVal, long double& recvVal,  MPI_Op op, int root)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 1, MPI_LONG_DOUBLE, op, root, getGlobalCommunicator());
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::reduce<__float128 >(__float128 sendVal, __float128& recvVal,  MPI_Op op, int root)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), sizeof(__float128), MPI_CHAR, op, root, getGlobalCommunicator());
}
#endif

template <>
void MpiManager::reduce<Complex<double> >(Complex<double> sendVal, Complex<double>& recvVal,  MPI_Op op, int root)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 2, MPI_DOUBLE, op, root, getGlobalCommunicator());
}

template <>
void MpiManager::reduce<Complex<long double> >(Complex<long double> sendVal, Complex<long double>& recvVal,  MPI_Op op, int root)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 2, MPI_LONG_DOUBLE, op, root, getGlobalCommunicator());
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::reduce<Complex<__float128> >(Complex<__float128> sendVal, Complex<__float128>& recvVal,  MPI_Op op, int root)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 2*sizeof(__float128), MPI_CHAR, op, root, getGlobalCommunicator());
}
#endif

template <>
void MpiManager::reduceVect<char>(std::vector<char>& sendVal,
                                  std::vector<char>& recvVal, MPI_Op op, int root)
{
    if (!ok) return;
    if (sendVal.empty()) return;
    std::vector<int> tmpSend(sendVal.begin(), sendVal.end());
    std::vector<int> tmpRecv;
    void* recvPtr;
    if (getRank()==root) {
        tmpRecv.resize(tmpSend.size());
        recvPtr = static_cast<void*>(&(tmpRecv[0]));
    }
    else {
        recvPtr = 0;
    }
    MPI_Reduce(static_cast<void*>(&(tmpSend[0])), recvPtr,
               sendVal.size(), MPI_INT, op, root, getGlobalCommunicator());
    recvVal.assign(tmpRecv.begin(), tmpRecv.end());
}

template <>
void MpiManager::reduceVect<int>(std::vector<int>& sendVal,
                                 std::vector<int>& recvVal, MPI_Op op, int root)
{
    if (!ok) return;
    if (sendVal.empty()) return;
    void* recvPtr;
    if (getRank()==root) {
        recvVal.resize(sendVal.size());
        recvPtr = static_cast<void*>(&(recvVal[0]));
    }
    else {
        recvPtr = 0;
    }
    MPI_Reduce(static_cast<void*>(&(sendVal[0])), recvPtr,
               sendVal.size(), MPI_INT, op, root, getGlobalCommunicator());
}

template <>
void MpiManager::reduceVect<long>(std::vector<long>& sendVal,
                                  std::vector<long>& recvVal, MPI_Op op, int root)
{
    if (!ok) return;
    if (sendVal.empty()) return;
    void* recvPtr;
    if (getRank()==root) {
        recvVal.resize(sendVal.size());
        recvPtr = static_cast<void*>(&(recvVal[0]));
    }
    else {
        recvPtr = 0;
    }
    MPI_Reduce(static_cast<void*>(&(sendVal[0])), recvPtr,
               sendVal.size(), MPI_LONG, op, root, getGlobalCommunicator());
}

template <>
void MpiManager::reduceVect<float>(std::vector<float>& sendVal,
                                   std::vector<float>& recvVal, MPI_Op op, int root)
{
    if (!ok) return;
    if (sendVal.empty()) return;
    void* recvPtr;
    if (getRank()==root) {
        recvVal.resize(sendVal.size());
        recvPtr = static_cast<void*>(&(recvVal[0]));
    }
    else {
        recvPtr = 0;
    }
    MPI_Reduce(static_cast<void*>(&(sendVal[0])), recvPtr,
               sendVal.size(), MPI_FLOAT, op, root, getGlobalCommunicator());
}

template <>
void MpiManager::reduceVect<double>(std::vector<double>& sendVal,
                                    std::vector<double>& recvVal, MPI_Op op, int root)
{
    if (!ok) return;
    if (sendVal.empty()) return;
    void* recvPtr;
    if (getRank()==root) {
        recvVal.resize(sendVal.size());
        recvPtr = static_cast<void*>(&(recvVal[0]));
    }
    else {
        recvPtr = 0;
    }
    MPI_Reduce(static_cast<void*>(&(sendVal[0])), recvPtr,
               sendVal.size(), MPI_DOUBLE, op, root, getGlobalCommunicator());
}

template <>
void MpiManager::reduceVect<long double>(std::vector<long double>& sendVal,
                                         std::vector<long double>& recvVal, MPI_Op op, int root)
{
    if (!ok) return;
    if (sendVal.empty()) return;
    void* recvPtr;
    if (getRank()==root) {
        recvVal.resize(sendVal.size());
        recvPtr = static_cast<void*>(&(recvVal[0]));
    }
    else {
        recvPtr = 0;
    }
    MPI_Reduce(static_cast<void*>(&(sendVal[0])), recvPtr,
               sendVal.size(), MPI_LONG_DOUBLE, op, root, getGlobalCommunicator());
}

#ifdef PLB_BGP
template <>
void MpiManager::reduceVect<long long>(std::vector<long long>& sendVal,
                                       std::vector<long long>& recvVal, MPI_Op op, int root)
{
    if (!ok) return;
    if (sendVal.empty()) return;
    void* recvPtr;
    if (getRank()==root) {
        recvVal.resize(sendVal.size());
        recvPtr = static_cast<void*>(&(recvVal[0]));
    }
    else {
        recvPtr = 0;
    }
    MPI_Reduce(static_cast<void*>(&(sendVal[0])), recvPtr,
               sendVal.size(), MPI_LONG_LONG, op, root, getGlobalCommunicator());
}
#endif

#ifdef PLB_USE_GCC
template <>
void MpiManager::reduceVect<__float128 >(std::vector<__float128 >& sendVal,
                                         std::vector<__float128 >& recvVal, MPI_Op op, int root)
{
    if (!ok) return;
    if (sendVal.empty()) return;
    void* recvPtr;
    if (getRank()==root) {
        recvVal.resize(sendVal.size());
        recvPtr = static_cast<void*>(&(recvVal[0]));
    }
    else {
        recvPtr = 0;
    }
    MPI_Reduce(static_cast<void*>(&(sendVal[0])), recvPtr,
               sendVal.size()*sizeof(__float128), MPI_CHAR, op, root, getGlobalCommunicator());
}
#endif

template <>
void MpiManager::reduceVect<Complex<double> >(std::vector<Complex<double> >& sendVal,
                                              std::vector<Complex<double> >& recvVal, MPI_Op op, int root)
{
    if (!ok) return;
    if (sendVal.empty()) return;
    void* recvPtr;
    if (getRank()==root) {
        recvVal.resize(sendVal.size());
        recvPtr = static_cast<void*>(&(recvVal[0]));
    }
    else {
        recvPtr = 0;
    }
    MPI_Reduce(static_cast<void*>(&(sendVal[0])), recvPtr,
               2*sendVal.size(), MPI_DOUBLE, op, root, getGlobalCommunicator());
}

template <>
void MpiManager::reduceVect<Complex<long double> >(std::vector<Complex<long double> >& sendVal,
                                                   std::vector<Complex<long double> >& recvVal, MPI_Op op, int root)
{
    if (!ok) return;
    if (sendVal.empty()) return;
    void* recvPtr;
    if (getRank()==root) {
        recvVal.resize(sendVal.size());
        recvPtr = static_cast<void*>(&(recvVal[0]));
    }
    else {
        recvPtr = 0;
    }
    MPI_Reduce(static_cast<void*>(&(sendVal[0])), recvPtr,
               2*sendVal.size(), MPI_LONG_DOUBLE, op, root, getGlobalCommunicator());
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::reduceVect<Complex<__float128> >(std::vector<Complex<__float128> >& sendVal,
                                                  std::vector<Complex<__float128> >& recvVal, MPI_Op op, int root)
{
    if (!ok) return;
    if (sendVal.empty()) return;
    void* recvPtr;
    if (getRank()==root) {
        recvVal.resize(sendVal.size());
        recvPtr = static_cast<void*>(&(recvVal[0]));
    }
    else {
        recvPtr = 0;
    }
    MPI_Reduce(static_cast<void*>(&(sendVal[0])), recvPtr,
               2*sendVal.size(), MPI_CHAR, op, root, getGlobalCommunicator());
}
#endif

template <>
void MpiManager::allReduceVect<char>(std::vector<char>& sendRecvVal, MPI_Op op)
{
    if (!ok) return;
    if (sendRecvVal.empty()) return;
    std::vector<char> recvVal(sendRecvVal.size());
    MPI_Allreduce( static_cast<void*>(&(sendRecvVal[0])),
                   static_cast<void*>(&(recvVal[0])),
                   sendRecvVal.size(), MPI_CHAR, op, getGlobalCommunicator() );
    sendRecvVal.swap(recvVal);
}

template <>
void MpiManager::allReduceVect<int>(std::vector<int>& sendRecvVal, MPI_Op op)
{
    if (!ok) return;
    if (sendRecvVal.empty()) return;
    std::vector<int> recvVal(sendRecvVal.size());
    MPI_Allreduce( static_cast<void*>(&(sendRecvVal[0])),
                   static_cast<void*>(&(recvVal[0])),
                   sendRecvVal.size(), MPI_INT, op, getGlobalCommunicator() );
    sendRecvVal.swap(recvVal);
}

template <>
void MpiManager::allReduceVect<long>(std::vector<long>& sendRecvVal, MPI_Op op)
{
    if (!ok) return;
    if (sendRecvVal.empty()) return;
    std::vector<long> recvVal(sendRecvVal.size());
    MPI_Allreduce( static_cast<void*>(&(sendRecvVal[0])),
                   static_cast<void*>(&(recvVal[0])),
                   sendRecvVal.size(), MPI_LONG, op, getGlobalCommunicator() );
    sendRecvVal.swap(recvVal);
}

#ifdef PLB_BGP
template <>
void MpiManager::allReduceVect<long long>(std::vector<long long>& sendRecvVal, MPI_Op op)
{
    if (!ok) return;
    if (sendRecvVal.empty()) return;
    std::vector<long long> recvVal(sendRecvVal.size());
    MPI_Allreduce( static_cast<void*>(&(sendRecvVal[0])),
                   static_cast<void*>(&(recvVal[0])),
                   sendRecvVal.size(), MPI_LONG_LONG, op, getGlobalCommunicator() );
    sendRecvVal.swap(recvVal);
}

template <>
void MpiManager::allReduceVect<unsigned long long>(std::vector<unsigned long long>& sendRecvVal, MPI_Op op)
{
    if (!ok) return;
    if (sendRecvVal.empty()) return;
    std::vector<unsigned long long> recvVal(sendRecvVal.size());
    MPI_Allreduce( static_cast<void*>(&(sendRecvVal[0])),
                   static_cast<void*>(&(recvVal[0])),
                   sendRecvVal.size(), MPI_UNSIGNED_LONG_LONG, op, getGlobalCommunicator() );
    sendRecvVal.swap(recvVal);
}

#endif

template <>
void MpiManager::allReduceVect<float>(std::vector<float>& sendRecvVal, MPI_Op op)
{
    if (!ok) return;
    if (sendRecvVal.empty()) return;
    std::vector<float> recvVal(sendRecvVal.size());
    MPI_Allreduce( static_cast<void*>(&(sendRecvVal[0])),
                   static_cast<void*>(&(recvVal[0])),
                   sendRecvVal.size(), MPI_FLOAT, op, getGlobalCommunicator() );
    sendRecvVal.swap(recvVal);
}

template <>
void MpiManager::allReduceVect<double>(std::vector<double>& sendRecvVal, MPI_Op op)
{
    if (!ok) return;
    if (sendRecvVal.empty()) return;
    std::vector<double> recvVal(sendRecvVal.size());
    MPI_Allreduce( static_cast<void*>(&(sendRecvVal[0])),
                   static_cast<void*>(&(recvVal[0])),
                   sendRecvVal.size(), MPI_DOUBLE, op, getGlobalCommunicator() );
    sendRecvVal.swap(recvVal);
}

template <>
void MpiManager::allReduceVect<long double>(std::vector<long double>& sendRecvVal, MPI_Op op)
{
    if (!ok) return;
    if (sendRecvVal.empty()) return;
    std::vector<long double> recvVal(sendRecvVal.size());
    MPI_Allreduce( static_cast<void*>(&(sendRecvVal[0])),
                   static_cast<void*>(&(recvVal[0])),
                   sendRecvVal.size(), MPI_LONG_DOUBLE, op, getGlobalCommunicator() );
    sendRecvVal.swap(recvVal);
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::allReduceVect<__float128 >(std::vector<__float128 >& sendRecvVal, MPI_Op op)
{
    // Reductions are not defined for this type.
    PLB_ASSERT( false );
}
#endif

template <>
void MpiManager::allReduceVect<Complex<double> >(std::vector<Complex<double> >& sendRecvVal, MPI_Op op)
{
    if (!ok) return;
    if (sendRecvVal.empty()) return;
    std::vector<Complex<double> > recvVal(sendRecvVal.size());
    MPI_Allreduce( static_cast<void*>(&(sendRecvVal[0])),
                   static_cast<void*>(&(recvVal[0])),
                   2*sendRecvVal.size(), MPI_DOUBLE, op, getGlobalCommunicator() );
    sendRecvVal.swap(recvVal);
}

template <>
void MpiManager::allReduceVect<Complex<long double> >(std::vector<Complex<long double> >& sendRecvVal, MPI_Op op)
{
    if (!ok) return;
    if (sendRecvVal.empty()) return;
    std::vector<Complex<long double> > recvVal(sendRecvVal.size());
    MPI_Allreduce( static_cast<void*>(&(sendRecvVal[0])),
                   static_cast<void*>(&(recvVal[0])),
                   2*sendRecvVal.size(), MPI_LONG_DOUBLE, op, getGlobalCommunicator() );
    sendRecvVal.swap(recvVal);
}

#ifdef PLB_USE_GCC
template <>
void MpiManager::allReduceVect<Complex<__float128> >(std::vector<Complex<__float128> >& sendRecvVal, MPI_Op op)
{
    // Reductions are not defined for this type.
    PLB_ASSERT( false );
}
#endif

template <>
void MpiManager::reduceAndBcast<char>(char& reductVal, MPI_Op op, int root)
{
    if (!ok) return;
    char recvVal;
    MPI_Reduce(&reductVal, &recvVal, 1, MPI_CHAR, op, root, getGlobalCommunicator());
    reductVal = recvVal;
    MPI_Bcast(&reductVal, 1, MPI_CHAR, root, getGlobalCommunicator());

}

template <>
void MpiManager::reduceAndBcast<int>(int& reductVal, MPI_Op op, int root)
{
    if (!ok) return;
    int recvVal;
    MPI_Reduce(&reductVal, &recvVal, 1, MPI_INT, op, root, getGlobalCommunicator());
    reductVal = recvVal;
    MPI_Bcast(&reductVal, 1, MPI_INT, root, getGlobalCommunicator());

}

template <>
void MpiManager::reduceAndBcast<long>(long& reductVal, MPI_Op op, int root)
{
    if (!ok) return;
    long recvVal;
    MPI_Reduce(&reductVal, &recvVal, 1, MPI_LONG, op, root, getGlobalCommunicator());
    reductVal = recvVal;
    MPI_Bcast(&reductVal, 1, MPI_LONG, root, getGlobalCommunicator());

}

#ifdef PLB_BGP
template <>
void MpiManager::reduceAndBcast<long long>(long long& reductVal, MPI_Op op, int root)
{
    if (!ok) return;
    int recvVal;
    MPI_Reduce(&reductVal, &recvVal, sizeof(long long), MPI_CHAR, op, root, getGlobalCommunicator());
    reductVal = recvVal;
    MPI_Bcast(&reductVal, sizeof(long long), MPI_CHAR, root, getGlobalCommunicator());

}

template <>
void MpiManager::reduceAndBcast<unsigned long long>(unsigned long long& reductVal, MPI_Op op, int root)
{
    if (!ok) return;
    int recvVal;
    MPI_Reduce(&reductVal, &recvVal, sizeof(unsigned long long), MPI_CHAR, op, root, getGlobalCommunicator());
    reductVal = recvVal;
    MPI_Bcast(&reductVal, sizeof(unsigned long long), MPI_CHAR, root, getGlobalCommunicator());

}
#endif

template <>
void MpiManager::reduceAndBcast<float>(float& reductVal, MPI_Op op, int root)
{
    if (!ok) return;
    float recvVal;
    MPI_Reduce(&reductVal, &recvVal, 1, MPI_FLOAT, op, root, getGlobalCommunicator());
    reductVal = recvVal;
    MPI_Bcast(&reductVal, 1, MPI_FLOAT, root, getGlobalCommunicator());

}

template <>
void MpiManager::reduceAndBcast<double>(double& reductVal, MPI_Op op, int root)
{
    if (!ok) return;
    double recvVal;
    MPI_Reduce(&reductVal, &recvVal, 1, MPI_DOUBLE, op, root, getGlobalCommunicator());
    reductVal = recvVal;
    MPI_Bcast(&reductVal, 1, MPI_DOUBLE, root, getGlobalCommunicator());

}

template <>
void MpiManager::reduceAndBcast<long double>(long double& reductVal, MPI_Op op, int root)
{
    if (!ok) return;
    long double recvVal;
    MPI_Reduce(&reductVal, &recvVal, 1, MPI_LONG_DOUBLE, op, root, getGlobalCommunicator());
    reductVal = recvVal;
    MPI_Bcast(&reductVal, 1, MPI_LONG_DOUBLE, root, getGlobalCommunicator());

}

#ifdef PLB_USE_GCC
template <>
void MpiManager::reduceAndBcast<__float128 >(__float128& reductVal, MPI_Op op, int root)
{
    if (!ok) return;
    __float128  recvVal;
    MPI_Reduce(&reductVal, &recvVal, sizeof(__float128), MPI_CHAR, op, root, getGlobalCommunicator());
    reductVal = recvVal;
    MPI_Bcast(&reductVal, sizeof(__float128), MPI_CHAR, root, getGlobalCommunicator());

}
#endif

template <>
void MpiManager::reduceAndBcast<Complex<double> >(Complex<double>& reductVal, MPI_Op op, int root)
{
    if (!ok) return;
    Complex<double>  recvVal;
    MPI_Reduce(&reductVal, &recvVal, 2, MPI_DOUBLE, op, root, getGlobalCommunicator());
    reductVal = recvVal;
    MPI_Bcast(&reductVal, 2, MPI_DOUBLE, root, getGlobalCommunicator());

}

template <>
void MpiManager::reduceAndBcast<Complex<long double> >(Complex<long double>& reductVal, MPI_Op op, int root)
{
    if (!ok) return;
    Complex<long double>  recvVal;
    MPI_Reduce(&reductVal, &recvVal, 2, MPI_LONG_DOUBLE, op, root, getGlobalCommunicator());
    reductVal = recvVal;
    MPI_Bcast(&reductVal, 2, MPI_DOUBLE, root, getGlobalCommunicator());

}

#ifdef PLB_USE_GCC
template <>
void MpiManager::reduceAndBcast<Complex<__float128> >(Complex<__float128>& reductVal, MPI_Op op, int root)
{
    if (!ok) return;
    Complex<__float128>  recvVal;
    MPI_Reduce(&reductVal, &recvVal, 2*sizeof(__float128), MPI_CHAR, op, root, getGlobalCommunicator());
    reductVal = recvVal;
    MPI_Bcast(&reductVal, 2*sizeof(__float128), MPI_CHAR, root, getGlobalCommunicator());

}
#endif

void MpiManager::wait(MPI_Request* request, MPI_Status* status)
{
    if (!ok) return;
    MPI_Wait(request, status);
}

}  // namespace global

}  // namespace plb

#endif  // PLB_MPI_PARALLEL
