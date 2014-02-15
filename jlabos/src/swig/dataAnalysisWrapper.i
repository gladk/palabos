namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * dataProcessors/dataAnalysisWrapper.i
 */

%{
#include "JLABOS_ROOT/plbWrapper/block/dataAnalysisWrapper2D.h"
#include "JLABOS_ROOT/plbWrapper/block/dataAnalysisWrapper3D.h"
#include "PALABOS_ROOT/src/dataProcessors/ntensorAnalysisWrapper2D.h"
#include "PALABOS_ROOT/src/dataProcessors/ntensorAnalysisWrapper3D.h"
%}


namespace plb {
    %include "arrays_java.i";
    %apply double[] {double *};
    %newobject copyConvert;
    %newobject maskedCopyConvert;
    %newobject extractComponent;
    %newobject maskedExtractComponent;
    %newobject computeNorm;
    %newobject maskedComputeNorm;
    %newobject computeNormSqr;
    %newobject maskedComputeNormSqr;
    %newobject computeSymmetricTensorNorm;
    %newobject maskedComputeSymmetricTensorNorm;
    %newobject computeSymmetricTensorNormSqr;
    %newobject maskedComputeSymmetricTensorNormSqr;
    %newobject computeSymmetricTensorTrace;
    %newobject maskedComputeSymmetricTensorTrace;
    %newobject computeVorticity;
    %newobject maskedComputeVorticity;
    %newobject computeBulkVorticity;
    %newobject maskedComputeBulkVorticity;
    %newobject computeStrainRate;
    %newobject maskedComputeStrainRate;
    %newobject computeStrainRate;
    %newobject maskedComputeStrainRate;
    %newobject add;
    %newobject maskedAdd;
    %newobject subtract;
    %newobject maskedSubtract;
    %newobject multiply;
    %newobject maskedMultiply;
    %newobject divide;
    %newobject maskedDivide;
    %newobject toThePower;
    %newobject maskedToThePower;
    %newobject equals;
    %newobject maskedEquals;
    %newobject lessThan;
    %newobject maskedLessThan;
    %newobject lessEqual;
    %newobject maskedLessEqual;
    %newobject greaterThan;
    %newobject maskedGreaterThan;
    %newobject greaterEqual;
    %newobject maskedGreaterEqual;
    %newobject logicalAnd;
    %newobject maskedLogicalAnd;
    %newobject logicalOr;
    %newobject maskedLogicalOr;
    %newobject negate;
    %newobject maskedNegate;
}

%include "JLABOS_ROOT/plbWrapper/block/dataAnalysisWrapper2D.h"
%include "JLABOS_ROOT/plbWrapper/block/dataAnalysisWrapper3D.h"
%include "PALABOS_ROOT/src/dataProcessors/ntensorAnalysisWrapper2D.h"
%include "PALABOS_ROOT/src/dataProcessors/ntensorAnalysisWrapper3D.h"


/*%apply (PRECOMP_T* INPLACE_ARRAY1, int DIM1) {(PRECOMP_T* result, int size)};*/
/*%apply (PRECOMP_T* IN_ARRAY1, int DIM1) {(PRECOMP_T* result, int size)};*/
%template(PRECOMP_T_plbAverage) plb::computeAverage<PRECOMP_T>;
%template(PRECOMP_T_plbBoundedAverage) plb::computeBoundedAverage<PRECOMP_T>;
%template(PRECOMP_T_plbMin) plb::computeMin<PRECOMP_T>;
%template(PRECOMP_T_plbMax) plb::computeMax<PRECOMP_T>;
%template(PRECOMP_T_plbCopyConvertInt) plb::copyConvert<PRECOMP_T,int>;
%template(PRECOMP_T_plbCopyConvertFloat) plb::copyConvert<PRECOMP_T,float>;
%template(PRECOMP_T_plbCopyConvertDouble) plb::copyConvert<PRECOMP_T,double>;
%template(PRECOMP_T_plbExtractComponent) plb::extractComponent<PRECOMP_T>;
%template(PRECOMP_T_plbComputeNorm) plb::computeNorm<PRECOMP_T>;
%template(PRECOMP_T_plbComputeNormSqr) plb::computeNormSqr<PRECOMP_T>;
%template(PRECOMP_T_plbComputeSymmetricTensorNorm) plb::computeSymmetricTensorNorm<PRECOMP_T>;
%template(PRECOMP_T_plbComputeSymmetricTensorNormSqr) plb::computeSymmetricTensorNormSqr<PRECOMP_T>;
%template(PRECOMP_T_plbComputeSymmetricTensorTrace) plb::computeSymmetricTensorTrace<PRECOMP_T>;
%template(PRECOMP_T_plbComputeVorticity) plb::computeVorticity<PRECOMP_T>;
%template(PRECOMP_T_plbComputeBulkVorticity) plb::computeBulkVorticity<PRECOMP_T>;
%template(PRECOMP_T_plbComputeStrainRate) plb::computeStrainRate<PRECOMP_T>;
%template(PRECOMP_T_plbComputeBulkStrainRate) plb::computeBulkStrainRate<PRECOMP_T>;
%template(PRECOMP_T_plbUPO_ScalarProduct) plb::compute_UPO_ScalarProduct<PRECOMP_T>;

%template(PRECOMP_T_m_plbAverage) plb::maskedComputeAverage<PRECOMP_T>;
%template(PRECOMP_T_m_plbBoundedAverage) plb::maskedComputeBoundedAverage<PRECOMP_T>;
%template(PRECOMP_T_m_plbMin) plb::maskedComputeMin<PRECOMP_T>;
%template(PRECOMP_T_m_plbMax) plb::maskedComputeMax<PRECOMP_T>;
%template(PRECOMP_T_m_plbCopyConvertInt) plb::maskedCopyConvert<PRECOMP_T,int>;
%template(PRECOMP_T_m_plbCopyConvertFloat) plb::maskedCopyConvert<PRECOMP_T,float>;
%template(PRECOMP_T_m_plbCopyConvertDouble) plb::maskedCopyConvert<PRECOMP_T,double>;
%template(PRECOMP_T_m_plbExtractComponent) plb::maskedExtractComponent<PRECOMP_T>;
%template(PRECOMP_T_m_plbComputeNorm) plb::maskedComputeNorm<PRECOMP_T>;
%template(PRECOMP_T_m_plbComputeNormSqr) plb::maskedComputeNormSqr<PRECOMP_T>;
%template(PRECOMP_T_m_plbComputeSymmetricTensorNorm) plb::maskedComputeSymmetricTensorNorm<PRECOMP_T>;
%template(PRECOMP_T_m_plbComputeSymmetricTensorNormSqr) plb::maskedComputeSymmetricTensorNormSqr<PRECOMP_T>;
%template(PRECOMP_T_m_plbComputeSymmetricTensorTrace) plb::maskedComputeSymmetricTensorTrace<PRECOMP_T>;
%template(PRECOMP_T_m_plbComputeVorticity) plb::maskedComputeVorticity<PRECOMP_T>;
%template(PRECOMP_T_m_plbComputeBulkVorticity) plb::maskedComputeBulkVorticity<PRECOMP_T>;
%template(PRECOMP_T_m_plbComputeStrainRate) plb::maskedComputeStrainRate<PRECOMP_T>;
%template(PRECOMP_T_m_plbComputeBulkStrainRate) plb::maskedComputeBulkStrainRate<PRECOMP_T>;
%template(PRECOMP_T_m_plbUPO_ScalarProduct) plb::masked_compute_UPO_ScalarProduct<PRECOMP_T>;

/*%apply (PRECOMP_T* IN_ARRAY1, int DIM1) {(PRECOMP_T* alpha, int size)};*/

%template(PRECOMP_T_iassign)      plb::assignInPlace<PRECOMP_T>;
%template(PRECOMP_T_iadd)         plb::addInPlace<PRECOMP_T>;
%template(PRECOMP_T_add)          plb::add<PRECOMP_T>;
%template(PRECOMP_T_isubtract)    plb::subtractInPlace<PRECOMP_T>;
%template(PRECOMP_T_subtract)     plb::subtract<PRECOMP_T>;
%template(PRECOMP_T_imultiply)    plb::multiplyInPlace<PRECOMP_T>;
%template(PRECOMP_T_multiply)     plb::multiply<PRECOMP_T>;
%template(PRECOMP_T_idivide)      plb::divideInPlace<PRECOMP_T>;
%template(PRECOMP_T_divide)       plb::divide<PRECOMP_T>;
%template(PRECOMP_T_itoThePower)  plb::toThePowerInPlace<PRECOMP_T>;
%template(PRECOMP_T_toThePower)   plb::toThePower<PRECOMP_T>;
%template(PRECOMP_T_equals)       plb::equals<PRECOMP_T>;
%template(PRECOMP_T_lessThan)     plb::lessThan<PRECOMP_T>;
%template(PRECOMP_T_lessEqual)    plb::lessEqual<PRECOMP_T>;
%template(PRECOMP_T_greaterThan)  plb::greaterThan<PRECOMP_T>;
%template(PRECOMP_T_greaterEqual) plb::greaterEqual<PRECOMP_T>;
%template(PRECOMP_T_logicalAnd)   plb::logicalAnd<PRECOMP_T>;
%template(PRECOMP_T_logicalOr)    plb::logicalOr<PRECOMP_T>;
%template(PRECOMP_T_negate)       plb::negate<PRECOMP_T>;

%template(PRECOMP_T_m_iassign)      plb::maskedAssignInPlace<PRECOMP_T>;
%template(PRECOMP_T_m_iadd)         plb::maskedAddInPlace<PRECOMP_T>;
%template(PRECOMP_T_m_add)          plb::maskedAdd<PRECOMP_T>;
%template(PRECOMP_T_m_isubtract)    plb::maskedSubtractInPlace<PRECOMP_T>;
%template(PRECOMP_T_m_subtract)     plb::maskedSubtract<PRECOMP_T>;
%template(PRECOMP_T_m_imultiply)    plb::maskedMultiplyInPlace<PRECOMP_T>;
%template(PRECOMP_T_m_multiply)     plb::maskedMultiply<PRECOMP_T>;
%template(PRECOMP_T_m_idivide)      plb::maskedDivideInPlace<PRECOMP_T>;
%template(PRECOMP_T_m_divide)       plb::maskedDivide<PRECOMP_T>;
%template(PRECOMP_T_m_itoThePower)  plb::maskedToThePowerInPlace<PRECOMP_T>;
%template(PRECOMP_T_m_toThePower)   plb::maskedToThePower<PRECOMP_T>;
%template(PRECOMP_T_m_equals)       plb::maskedEquals<PRECOMP_T>;
%template(PRECOMP_T_m_lessThan)     plb::maskedLessThan<PRECOMP_T>;
%template(PRECOMP_T_m_lessEqual)    plb::maskedLessEqual<PRECOMP_T>;
%template(PRECOMP_T_m_greaterThan)  plb::maskedGreaterThan<PRECOMP_T>;
%template(PRECOMP_T_m_greaterEqual) plb::maskedGreaterEqual<PRECOMP_T>;
%template(PRECOMP_T_m_logicalAnd)   plb::maskedLogicalAnd<PRECOMP_T>;
%template(PRECOMP_T_m_logicalOr)    plb::maskedLogicalOr<PRECOMP_T>;
%template(PRECOMP_T_m_negate)       plb::maskedNegate<PRECOMP_T>;

