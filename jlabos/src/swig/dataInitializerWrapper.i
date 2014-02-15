namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * dataProcessors/dataInitializerWrapper.i
 */

%{
#include "JLABOS_ROOT/plbWrapper/block/dataInitializerWrapper2D.h"
#include "JLABOS_ROOT/plbWrapper/block/dataInitializerWrapper3D.h"
%}

%include "JLABOS_ROOT/plbWrapper/block/dataInitializerWrapper2D.h"
%include "JLABOS_ROOT/plbWrapper/block/dataInitializerWrapper3D.h"

/*%include array.i*/
%apply int[] {int *}
%template(PRECOMP_T_plbSetToConstant)    plb::setToConstant<PRECOMP_T>;
%template(PRECOMP_T_plbSetToCoordinate)  plb::setToCoordinate<PRECOMP_T>;
%template(PRECOMP_T_plbSetToCoordinates) plb::setToCoordinates<PRECOMP_T>;
%template(PRECOMP_T_plbAssignComponent)  plb::assignComponent<PRECOMP_T>;

%template(PRECOMP_T_m_plbSetToConstant)    plb::maskedSetToConstant<PRECOMP_T>;
%template(PRECOMP_T_m_plbSetToCoordinate)  plb::maskedSetToCoordinate<PRECOMP_T>;
%template(PRECOMP_T_m_plbSetToCoordinates) plb::maskedSetToCoordinates<PRECOMP_T>;
%template(PRECOMP_T_m_plbAssignComponent)  plb::maskedAssignComponent<PRECOMP_T>;

