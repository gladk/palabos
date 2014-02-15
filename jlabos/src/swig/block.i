namespace plb {
    template<typename T> class MultiNTensorField3D;
}

/*
 * block.i
 */

%module MODULE_NAME;

%{
#define SWIG_FILE_WITH_INIT
%}

namespace std {
    class exception;
}



namespace plb {
    typedef long plint;
    typedef size_t pluint;
    namespace IndexOrdering {
        enum OrderingT {forward, backward, memorySaving};
    }
}



/*%typemap(javaimports) Box2D "plb.core.Box2D";
%include "JLABOS_ROOT/swig/core/pre_processed/geometry2D.i";
%typemap(javaimports) Box3D "plb.core.Box3D";
%include "/home/ysagon/strep/yann/jlabos/src/swig/core/pre_processed/geometry3D.i";*/

namespace plb {
    %include "arrays_java.i";
    %apply double[] {double *};
    class PlbException;
    class PlbMemoryException;
    class PlbIOException;
    class PlbNetworkException;
    class PlbLogicException;
    class PlbOutOfRangeException;

    struct Box2D;
    struct Box3D;

}

/* Rename namespace global, because "global" is
 * a protected keyword in Python.
 */
%rename plb::global plb_global;
%include "core_geometry2D.i";
%include "core_geometry3D.i";
%include "PRECOMP_T_block_multiDataField.i";
%include "PRECOMP_T_block_multiBlockGenerator.i";
%include "PRECOMP_T_block_dataInitializerWrapper.i";
%include "PRECOMP_T_block_dataAnalysisWrapper.i";

