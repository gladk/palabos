namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * blockLattice3d/dynamics.i
 */

%{
#include "PALABOS_ROOT/src/latticeBoltzmann/nearestNeighborLattices3D.h"
#include "PALABOS_ROOT/src/latticeBoltzmann/nearestNeighborLattices3D.hh"
#include "PALABOS_ROOT/src/core/dynamics.h"
%}


namespace plb {

template<typename T, class Descriptor>
struct Dynamics {
    virtual ~Dynamics();
    virtual int getId() const =0;
};

template<typename T, class Descriptor>
class BounceBack : public Dynamics<T,Descriptor> {
public:
    BounceBack(T virtualDensity);
    virtual int getId() const;
};

}

%template(FLOAT_T_DESCRIPTOR_3D_PlbDynamics) plb::Dynamics<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_PlbBounceBack) plb::BounceBack<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;

