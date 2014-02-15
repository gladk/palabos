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
# Leading directory of the SWIG source code
jlabosRoot=$(PALABOS_ROOT)/jlabos/src/

# main class name without extension
mainClass=Cavity2d
# Set optimization flags on/off
#
optimize     = true
# Set debug mode and debug flags on/off
debug        = true
# Set profiling flags on/off
profile      = false
# Java compiler
javaCompiler = javac

# Additional include paths
#includePaths =

##########################################################################
# All code below this line is just about forwarding the options
# to SConstruct. It is recommended not to modify anything there.
##########################################################################

SCons     = $(palabosRoot)/scons/scons.py -f ./SConstruct

SConsArgs = palabosRoot=$(palabosRoot) \
            jlabosRoot=$(jlabosRoot) \
            mainClass=$(mainClass) \
            optimize=$(optimize) \
            debug=$(debug) \
            profile=$(profile) \
            javaCompiler=$(javaCompiler) 

run: compilation
	mpirun -np 4 java -ea -jar -Djava.library.path=../../src/jlabos/plb:/usr/lib $(mainClass).jar
	
compilation:
	python $(SCons) $(SConsArgs)

clean:
	python $(SCons) -c $(SConsArgs)
