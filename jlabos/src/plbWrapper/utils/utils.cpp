/* utils for palabos wrapper */

#include <dlfcn.h>
#include <mpi.h>
#include <iostream>
#include "utils.h"
char* Mpi::loadMpi(void){
        MPISIZE = 0; 
        MPIRANK = 0;
        //int MPISIZE = 0; 
        //char **MPIRANK = 0;
        //dlopen("libmpi.so.0", RTLD_LAZY);
        dlopen("libmpi.so.0", RTLD_NOW | RTLD_GLOBAL);
        MPI_Init(&MPISIZE, &MPIRANK);
        //int nprocs;
        //MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
        //std::cout << "nb procs :" << nprocs << std::endl;
        return dlerror();
}

int Mpi::finalizeMpi(void){
	dlopen("libmpi.so.0", RTLD_NOW | RTLD_GLOBAL);
	return MPI_Finalize();
}

int Mpi::getRankMpi(void){
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	return rank;
}
