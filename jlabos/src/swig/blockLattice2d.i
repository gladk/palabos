namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * blockLattice2d.i
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

    /*%include "JLABOS_ROOT/swig/core/pre_processed/globalDefs.i"*/
    /*typedef ptrdiff_t plint;
    typedef size_t pluint;*/

}
/*
%typemap(javaimports) Box2D "plb.core.Box2D"
%include "JLABOS_ROOT/swig/core/pre_processed/geometry2D.i";
%typemap(javaimports) Box3D "plb.core.Box3D"
%include "JLABOS_ROOT/swig/core/pre_processed/geometry3D.i";
%typemap(javaimports) MultiNTensorField2D "plb.block.mydouble.MultiNTensorField2D"
%include "JLABOS_ROOT/swig/block/mydouble/multiDataField.i";
%typemap(javaimports) int_block "plb.block.myint"
*/

namespace plb {
    template<typename T> class MultiNTensorField2D;
}




namespace plb {
    class PlbException;
    class PlbMemoryException;
    class PlbIOException;
    class PlbNetworkException;
    class PlbLogicException;
    class PlbOutOfRangeException;

    struct Box2D;
}

/* Rename namespace global, because "global" is
 * a protected keyword in Python.
 */
%rename plb::global plb_global;

%include "FLOAT_T_d2q9_dynamics2D.i";
%include "FLOAT_T_d2q9_isoThermalDynamics2D.i";
%include "FLOAT_T_d2q9_multiBlockLattice2D.i";
%include "FLOAT_T_d2q9_multiBlockGenerator2D.i";
%include "FLOAT_T_d2q9_latticeInitializerWrapper2D.i";
%include "FLOAT_T_d2q9_latticeAnalysisWrapper2D.i";
%include "FLOAT_T_d2q9_boundaryCondition2D.i";
%include "FLOAT_T_d2q9_javaInterface2D.i";
%include "FLOAT_T_d2q9_javaMultiTensorFieldInterface2D.i";
%include "JLABOS_ROOT/swig/pre_processed/FLOAT_T_block_block.i";
