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

#include "core/plbDebug.h"
#include "parallelism/mpiManager.h"

#include <cerrno>
#include <cstdio>

#ifdef PLB_USE_POSIX
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

namespace plb {

void enableCoreDumps()
{
#if defined PLB_USE_POSIX && defined PLB_DEBUG
    struct rlimit lim;

    lim.rlim_cur = RLIM_INFINITY;
    lim.rlim_max = RLIM_INFINITY;

    if (setrlimit(RLIMIT_CORE, &lim)) {
        char buf[128];
        int myRank = (int) global::mpi().getRank();
        sprintf(buf, "enableCoreDumps(): setrlimit at process %d", myRank);
        perror(buf);
    }
#endif
}

void unbufferOutputStdStreams()
{
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);
}

}  // namespace plb
