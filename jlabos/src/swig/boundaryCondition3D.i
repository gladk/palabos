namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * plbWrapper/lattice/bcWrapper3D.i
 */

%{
#include "JLABOS_ROOT/plbWrapper/lattice/bcWrapper3D.h"
%}

namespace plb {

template<typename T, class Descriptor>
class OuterBoxBC {
public:
    void setVelocityCondition( MultiBlockLattice3D<T,Descriptor>& lattice,
                               Box3D domain );
    void setPressureCondition( MultiBlockLattice3D<T,Descriptor>& lattice,
                               Box3D domain );
private:
    OuterBoxBC();  /** Don't allow access to constructor from
                     * Python.
                     */
};

%newobject generateRegularizedBC;
template<typename T, class Descriptor>
OuterBoxBC<T,Descriptor>* generateRegularizedBC();

}


%template(FLOAT_T_DESCRIPTOR_3D_PlbOuterBoxBC) plb::OuterBoxBC<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_regularizedBC) plb::generateRegularizedBC<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
