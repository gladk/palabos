/*
 * core.i
 */

%module core;
%include "std_except.i";
%include "std_string.i";
%include "std_vector.i";

%exception %{
    try {
      $action
    } catch (std::bad_alloc &) {
      return $null;
    } catch (std::exception &e) {
      jclass clazz = jenv->FindClass("java/lang/Exception");
      jenv->ThrowNew(clazz, e.what());
      return $null;
    } catch (...) {
      jclass clazz = jenv->FindClass("java/lang/Exception");
      jenv->ThrowNew(clazz, "Unknown exception");
      return $null;
    }
%} 



typedef long plint;
typedef size_t pluint;


%template(StrVector) std::vector<std::string>;
%template(DoubleVector) std::vector<double>; 

/* Rename namespace global, because "global" is
 * a protected keyword in Python.
 */
%rename plb::global plb_global;

%include "core_globalDefs.i";
%include "core_runTimeDiagnostics.i";
%template(PlintVector) std::vector<plb::plint>;
%template(PluintVector) std::vector<plb::pluint>;
%include "core_plbInit.i";
%include "core_plbTimer.i";
%include "core_geometry2D.i"
%include "core_geometry3D.i"
%include "core_blockStatistics.i"

%{
#include "JLABOS_ROOT/plbWrapper/utils/utils.h"
%}
%include "JLABOS_ROOT/plbWrapper/utils/utils.h";

/*char* loadMpi(void);
int finalizeMpi(void);
bool isMPIRankZero();*/
