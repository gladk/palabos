namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * plbWrapper/lattice/bcWrapper2D.i
 */

%{
#include "JLABOS_ROOT/plbWrapper/lattice/bcWrapper2D.h"
%}

namespace plb {

template<typename T, class Descriptor>
class OuterBoxBC {
public:
    void setVelocityCondition( MultiBlockLattice2D<T,Descriptor>& lattice,
                               Box2D domain );
    void setPressureCondition( MultiBlockLattice2D<T,Descriptor>& lattice,
                               Box2D domain );
private:
    OuterBoxBC();  /** Don't allow access to constructor from
                     * Python.
                     */
};

%newobject generateRegularizedBC;
template<typename T, class Descriptor>
OuterBoxBC<T,Descriptor>* generateRegularizedBC();

}


%template(FLOAT_T_DESCRIPTOR_2D_PlbOuterBoxBC) plb::OuterBoxBC<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_regularizedBC) plb::generateRegularizedBC<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
