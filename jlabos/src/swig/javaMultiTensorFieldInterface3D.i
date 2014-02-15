namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * blockLattice2d/numPyInterface2D.i
 */

%{
#include "JLABOS_ROOT/plbWrapper/block/numPyInterface3D.h"
%}


namespace plb {

template<typename T> class MultiNTensorField3D;

%include "arrays_java.i";
%apply double[] {double *};
template<typename T>
class NTensorField2NumPy3D {
public:
    NTensorField2NumPy3D(MultiNTensorField3D<T>& field_);
    NTensorField2NumPy3D(MultiNTensorField3D<T>& field_, Box3D const& domain);
    void execute(T* array, int size);
    int getSize() const;
};
    

}  // namespace plb

%template(PRECOMP_T_NTensorFieldSerializer3D) plb::NTensorField2NumPy3D<PRECOMP_T>;
