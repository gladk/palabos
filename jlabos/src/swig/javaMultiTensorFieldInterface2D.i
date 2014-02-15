namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * blockLattice2d/numPyInterface2D.i
 */

%{
#include "JLABOS_ROOT/plbWrapper/block/numPyInterface2D.h"
%}


namespace plb {
%include "arrays_java.i";
%apply double[] {double *};
template<typename T>
class NTensorField2NumPy2D {
public:
    NTensorField2NumPy2D(MultiNTensorField2D<T>& field_);
    NTensorField2NumPy2D(MultiNTensorField2D<T>& field_, Box2D const& domain);
    void execute(T* array, int size);
    int getSize() const;
};

/* to be checked */
template<typename T>
class NumPy2NTensorField2D {
public:
    NumPy2NTensorField2D(MultiNTensorField2D<T>& field_);
    NumPy2NTensorField2D(MultiNTensorField2D<T>& field_, Box2D const& domain);
    void execute(T* array, int size);
    int getSize() const;
};


}  // namespace plb

/* to be checked */
%template(PRECOMP_T_NTensorFieldUnSerializer2D) plb::NumPy2NTensorField2D<PRECOMP_T>;


%template(PRECOMP_T_NTensorFieldSerializer2D) plb::NTensorField2NumPy2D<PRECOMP_T>;
