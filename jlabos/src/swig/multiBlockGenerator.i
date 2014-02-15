namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
  plbWrapper/multiBlockGenerators.i
 */

%{
#include "PALABOS_ROOT/src/multiBlock/multiBlockGenerator2D.h"
#include "PALABOS_ROOT/src/multiBlock/multiBlockGenerator3D.h"
%}

namespace plb {

%newobject generateMultiNTensorField2D;
template<typename T>
MultiNTensorField2D<T>* generateMultiNTensorField2D(Box2D const& domain, plint ndim);

%newobject generateNTensorFieldFromNTensor2D;
template<typename T1, typename T2>
MultiNTensorField2D<T2>*
    generateNTensorFieldFromNTensor2D (
    MultiNTensorField2D<T1> const& field,
    Box2D const& intersection, plint nDim );

%newobject generateMultiNTensorField3D;
template<typename T>
MultiNTensorField3D<T>* generateMultiNTensorField3D(Box3D const& domain, plint ndim);

%newobject generateNTensorFieldFromNTensor3D;
template<typename T1, typename T2>
MultiNTensorField3D<T2>*
    generateNTensorFieldFromNTensor3D (
    MultiNTensorField3D<T1> const& field,
    Box3D const& intersection, plint nDim );

}  // namespace plb

%template(PRECOMP_T_generateMultiNTensorField2D) plb::generateMultiNTensorField2D<PRECOMP_T>;
%template(PRECOMP_T_generateIntNTensorField2D)   plb::generateNTensorFieldFromNTensor2D<PRECOMP_T,int>;

%template(PRECOMP_T_generateMultiNTensorField3D) plb::generateMultiNTensorField3D<PRECOMP_T>;
%template(PRECOMP_T_generateIntNTensorField3D)   plb::generateNTensorFieldFromNTensor3D<PRECOMP_T,int>;

