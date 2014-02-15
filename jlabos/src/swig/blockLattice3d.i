namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * blockLattice3d.i
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
   /* typedef ptrdiff_t plint;
    typedef size_t pluint;*/

}
/*
%typemap(javaimports) Box2D "plb.core.Box2D"
%include "JLABOS_ROOT/swig/core/pre_processed/geometry2D.i";
%typemap(javaimports) Box3D "plb.core.Box3D"
%include "JLABOS_ROOT/swig/core/pre_processed/geometry3D.i";
%typemap(javaimports) MultiNTensorField3D "plb.block.mydouble.MultiNTensorField3D"
%include "JLABOS_ROOT/swig/block/mydouble/multiDataField.i";
*/
namespace plb {
    template<typename T> class MultiNTensorField3D;
}


namespace plb {
    class PlbException;
    class PlbMemoryException;
    class PlbIOException;
    class PlbNetworkException;
    class PlbLogicException;
    class PlbOutOfRangeException;

    struct Box3D;
}

/* Rename namespace global, because "global" is
 * a protected keyword in Python.
 */
%rename plb::global plb_global;


%include "core_geometry3D.i";
%include "core_geometry2D.i";
%include "FLOAT_T_block_multiDataField.i";
%include "MODULE_NAME_dynamics3D.i";
%include "MODULE_NAME_isoThermalDynamics3D.i";
%include "MODULE_NAME_multiBlockLattice3D.i";
%include "MODULE_NAME_multiBlockGenerator3D.i";
%include "MODULE_NAME_latticeInitializerWrapper3D.i";
%include "MODULE_NAME_latticeAnalysisWrapper3D.i";
%include "MODULE_NAME_boundaryCondition3D.i";
%include "MODULE_NAME_javaInterface3D.i";
%include "MODULE_NAME_javaMultiTensorFieldInterface3D.i";
