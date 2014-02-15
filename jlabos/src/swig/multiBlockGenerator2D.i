namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
  blockLattice2d/multiBlockGenerators2D.i
 */

%{
#include "PALABOS_ROOT/src/latticeBoltzmann/nearestNeighborLattices2D.h"
#include "PALABOS_ROOT/src/latticeBoltzmann/nearestNeighborLattices2D.hh"
#include "JLABOS_ROOT/plbWrapper/lattice/multiBlockGenerator2D.h"
%}

namespace plb {

%newobject generateMultiBlockLattice2D;
template<typename T, class Descriptor>
MultiBlockLattice2D<T,Descriptor>*
    generateMultiBlockLattice2D(Box2D const& domain, Dynamics<T,Descriptor> const* backgroundDynamics);

%newobject generateNTensorFieldFromLattice2D;
template<typename T1, class Descriptor, typename T2>
MultiNTensorField2D<T2>*
    generateNTensorFieldFromLattice2D (
    MultiBlockLattice2D<T1,Descriptor> const& lattice,
    Box2D const& domain, plint ndim );

}  // namespace plb

%template(FLOAT_T_DESCRIPTOR_2D_generateMultiBlockLattice) plb::generateMultiBlockLattice2D<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_generateIntNTensorField2D) plb::generateNTensorFieldFromLattice2D<FLOAT_T,plb::descriptors::DESCRIPTOR_2D,int>;
