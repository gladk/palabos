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

#include "core/globalDefs.h"
#include "parallelism/mpiManager.h"
#include "io/parallelIO.h"

namespace plb {

Parallel_referring_ostream pcout(std::cout);
Parallel_referring_ostream pcerr(std::cerr);
Parallel_referring_ostream pclog(std::clog);


/* *************** Class plb_ofstream ******************************** */

plb_ofstream::plb_ofstream()
    : devNullStream(&devNullBuffer),
      original (
          global::mpi().isMainProcessor() ?
            new std::ofstream : 0 )
{ } 

plb_ofstream::plb_ofstream(const char* filename, std::ostream::openmode mode)
    : devNullStream(&devNullBuffer),
      original (
          global::mpi().isMainProcessor() ?
            new std::ofstream(filename,mode) : 0 )
{ }

plb_ofstream::plb_ofstream(plb_ofstream const& rhs)
    : devNullStream(&devNullBuffer),
      original(0)
{ }

plb_ofstream& plb_ofstream::operator=(plb_ofstream const& rhs) {
    return *this;
}



plb_ofstream::~plb_ofstream() {
    delete original;
}

std::ostream& plb_ofstream::getOriginalStream()
{
    if (global::mpi().isMainProcessor()) {
        return *original;
    }
    else {
        return devNullStream;
    }
}

bool plb_ofstream::is_open() {
#ifdef PLB_MPI_PARALLEL
    int open = false;
    if (global::mpi().isMainProcessor()) {
        open = original->is_open();
    }
    global::mpi().bCast(&open, 1);
    return open;
#else
    return original->is_open();
#endif
}

void plb_ofstream::open(const char* filename, std::ostream::openmode mode)
{
    if (global::mpi().isMainProcessor()) {
        original->open(filename, mode);
    }
}

void plb_ofstream::close() {
    if (global::mpi().isMainProcessor()) {
        original->close();
    }
}


/* *************** Class plb_ifstream ******************************** */

plb_ifstream::plb_ifstream()
    : devNullStream(&devNullBuffer),
      original (
          global::mpi().isMainProcessor() ?
            new std::ifstream : 0 )
{ }

plb_ifstream::plb_ifstream(const char * filename, std::istream::openmode mode)
    : devNullStream(&devNullBuffer),
      original (
          global::mpi().isMainProcessor() ?
            new std::ifstream(filename,mode) : 0 )
{ }

plb_ifstream::plb_ifstream(plb_ifstream const& rhs)
    : devNullStream(&devNullBuffer),
      original(0)
{ }

plb_ifstream& plb_ifstream::operator=(plb_ifstream const& rhs) {
    return *this;
}

plb_ifstream::~plb_ifstream() {
    delete original;
}

std::istream& plb_ifstream::getOriginalStream() {
    if (global::mpi().isMainProcessor()) {
        return *original;
    }
    else {
        return devNullStream;
    }
}

bool plb_ifstream::is_open() {
#ifdef PLB_MPI_PARALLEL
    int open = false;
    if (global::mpi().isMainProcessor()) {
        open = original->is_open();
    }
    global::mpi().bCast(&open, 1);
    return open;
#else
    return original->is_open();
#endif
}

bool plb_ifstream::good() {
#ifdef PLB_MPI_PARALLEL
    int open = false;
    if (global::mpi().isMainProcessor()) {
        open = original->good();
    }
    global::mpi().bCast(&open, 1);
    return open;
#else
    return original->good();
#endif
}

void plb_ifstream::open(const char* filename, std::istream::openmode mode)
{
    if (global::mpi().isMainProcessor()) {
        original->open(filename, mode);
    }
}

void plb_ifstream::close() {
    if (global::mpi().isMainProcessor()) {
        original->close();
    }
}


}  // namespace plb
