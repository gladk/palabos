##########################################################################
## Makefile for the Jlabos C++ and SWIG based library.
##
## The present Makefile is a pure configuration file, in which 
## you can select compilation options. Compilation dependencies
## are managed automatically through the Python library SConstruct.
##
## If you don't have Python, or if compilation doesn't work for other
## reasons, consult the Palabos user's guide for instructions on manual
## compilation.
##########################################################################

# USE: multiple arguments are separated by spaces.
#   For example: projectFiles = file1.cpp file2.cpp
#                optimFlags   = -O -finline-functions

# Leading directory of the PALABOS source code
palabosRoot=$(PALABOS_ROOT)
# Leading directory of the SWIG source code (no slash at end)
jlabosRoot=$(PALABOS_ROOT)/jlabos/src

# Choose between generic and precompiled mode
precompiled  = true

# number of core used for compilation
nbOfCore     = 4
# Set optimization flags on/off
optimize     = true
# Set debug mode and debug flags on/off
debug        = true
# Set profiling flags on/off
profile      = false
# Set MPI-parallel mode on/off (parallelism in cluster-like environment)
MPIparallel  = true

# Set SMP-parallel mode on/off (shared-memory parallelism)
SMPparallel  = false
# Decide whether to include calls to the POSIX API. On non-POSIX systems,
#   including Windows, this flag must be false, unless a POSIX environment is
#   emulated (such as with Cygwin).
usePOSIX     = true
# Generate dynamic libraries.
dynamicLibrary = true

# Path to external libraries (other than Palabos)
libraryPaths =
# Swig executable
swigCompiler = swig
# Compiler to use without MPI parallelism
serialCXX    = g++
# Compiler to use with MPI parallelism
parallelCXX  = mpicxx
# General compiler flags (e.g. -Wall to turn on all warnings on g++)
compileFlags = -fno-strict-aliasing
# Compiler flags to use when optimization mode is on
optimFlags   = -O3
# Compiler flags to use when debug mode is on
debugFlags   = -g
# Compiler flags to use when profile mode is on
profileFlags = -pg
# Additional include paths
includePaths = /usr/lib/jvm/java-7-openjdk-amd64/include/ /usr/lib/jvm/java-7-openjdk-amd64/include/linux
# Additional linker flags
linkFlags = -Wl,--no-as-needed

##########################################################################
# All code below this line is just about forwarding the options
# to SConstruct. It is recommended not to modify anything there.
##########################################################################

SCons     = $(palabosRoot)/scons/scons.py -j $(nbOfCore) -f $(jlabosRoot)/SConstruct

SConsArgs = palabosRoot=$(palabosRoot) \
            jlabosRoot=$(jlabosRoot) \
            optimize=$(optimize) \
            debug=$(debug) \
	    swigCompiler=$(swigCompiler) \
            profile=$(profile) \
            MPIparallel=$(MPIparallel) \
            serialCXX=$(serialCXX) \
            parallelCXX=$(parallelCXX) \
            compileFlags="$(compileFlags)" \
            optimFlags="$(optimFlags)" \
            debugFlags="$(debugFlags)" \
	    profileFlags="$(profileFlags)" \
	    includePaths="$(includePaths)" \
	    linkFlags="$(linkFlags)"


SConsPalabos     = $(palabosRoot)/scons/scons.py -j $(nbOfCore) -f $(palabosRoot)/SConstruct

SConsArgsPalabos = palabosRoot=$(palabosRoot) \
            projectFiles="dummyMain.cpp" \
            precompiled=$(precompiled) \
            optimize=$(optimize) \
            debug=$(debug) \
            profile=$(profile) \
            MPIparallel=$(MPIparallel) \
            SMPparallel=$(SMPparallel) \
            usePOSIX=$(usePOSIX) \
            serialCXX=$(serialCXX) \
            parallelCXX=$(parallelCXX) \
            dynamicLibrary=$(dynamicLibrary) \
            compileFlags="$(compileFlags)" \
            linkFlags="$(linkFlags)" \
            optimFlags="$(optimFlags)" \
            debugFlags="$(debugFlags)" \
       profileFlags="$(profileFlags)" \
       libraryPaths="$(libraryPaths)" \
       includePaths="$(includePaths)" \
       libraries="$(libraries)"


compilation:
	cd compilePalabos; python $(SConsPalabos) $(SConsArgsPalabos)
	bash ./preprocess $(jlabosRoot)
	python $(SCons) $(SConsArgs)
	bash ./postprocess $(jlabosRoot) $(parallelCXX) $(linkFlags)

clean:
	python $(SCons) -c $(SConsArgs)
	cd compilePalabos; python $(SConsPalabos) -c $(SConsArgsPalabos)
	/bin/rm -f swig/tmp/*
	/bin/rm -f `find . -name "*~"`
	/bin/rm -f jlabos/plb/jlabos/*
	/bin/rm -f jlabos/plb/*.os
	/bin/rm -f swig/pre_processed/*
	/bin/rm -f precompiled/floatOnly/float/*
	/bin/rm -f precompiled/floatOnly/double/*
	/bin/rm -f precompiled/intOrFloat/int/*
	/bin/rm -f precompiled/intOrFloat/float/*
	/bin/rm -f precompiled/intOrFloat/double/*
	/bin/rm -f precompiled/lattice/d2q9_float/*
	/bin/rm -f precompiled/lattice/d2q9_double/*
	/bin/rm -f precompiled/lattice/d3q19_float/*
	/bin/rm -f precompiled/lattice/d3q19_double/*
	/bin/rm -f plbWrapper/block/int/*
	/bin/rm -f plbWrapper/block/float/*
	/bin/rm -f plbWrapper/block/double/*
	/bin/rm -f plbWrapper/lattice/d2q9_float/*
	/bin/rm -f plbWrapper/lattice/d2q9_double/*
	/bin/rm -f plbWrapper/lattice/d3q19_float/*
	/bin/rm -f plbWrapper/lattice/d3q19_double/*
