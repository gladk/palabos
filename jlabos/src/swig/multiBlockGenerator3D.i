namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
  blockLattice3d/multiBlockGenerators3D.i
 */

%{
#include "PALABOS_ROOT/src/latticeBoltzmann/nearestNeighborLattices3D.h"
#include "PALABOS_ROOT/src/latticeBoltzmann/nearestNeighborLattices3D.hh"
#include "JLABOS_ROOT/plbWrapper/lattice/multiBlockGenerator3D.h"
%}

namespace plb {

%newobject generateMultiBlockLattice3D;
template<typename T, class Descriptor>
MultiBlockLattice3D<T,Descriptor>*
    generateMultiBlockLattice3D(Box3D const& domain, Dynamics<T,Descriptor> const* backgroundDynamics);

%newobject generateNTensorFieldFromLattice3D;
template<typename T1, class Descriptor, typename T2>
MultiNTensorField3D<T2>*
    generateNTensorFieldFromLattice3D (
    MultiBlockLattice3D<T1,Descriptor> const& lattice,
    Box3D const& domain, plint ndim );

}  // namespace plb

%template(FLOAT_T_DESCRIPTOR_3D_generateMultiBlockLattice) plb::generateMultiBlockLattice3D<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_generateIntNTensorField3D)   plb::generateNTensorFieldFromLattice3D<FLOAT_T,plb::descriptors::DESCRIPTOR_3D,int>;
