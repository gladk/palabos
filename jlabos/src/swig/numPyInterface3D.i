namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * blockLattice3d/numPyInterface3D.i
 */

%{
#include "PALABOS_ROOT/src/latticeBoltzmann/nearestNeighborLattices3D.h"
#include "PALABOS_ROOT/src/latticeBoltzmann/nearestNeighborLattices3D.hh"
#include "JLABOS_ROOT/plbWrapper/lattice/numPyInterface3D.h"
%}

namespace plb {

template<typename T, class Descriptor>
class Lattice2NumPy3D {
public:
    Lattice2NumPy3D(MultiBlockLattice3D<T,Descriptor>& lattice_);
    Lattice2NumPy3D(MultiBlockLattice3D<T,Descriptor>& lattice_, Box3D const& domain);
    void execute(T* array, int size);
    int getSize() const;
};

template<typename T, class Descriptor>
class NumPy2Lattice3D {
public:
    NumPy2Lattice3D(MultiBlockLattice3D<T,Descriptor>& lattice_);
    NumPy2Lattice3D(MultiBlockLattice3D<T,Descriptor>& lattice_, Box3D const& domain);
    void execute(T* array, int size);
    int getSize() const;
};

}  // namespace plb

%apply (FLOAT_T* INPLACE_ARRAY1, int DIM1) {(FLOAT_T* array, int size)};
%template(FLOAT_T_DESCRIPTOR_3D_LatticeSerializer) plb::Lattice2NumPy3D<FLOAT_T, plb::descriptors::DESCRIPTOR_3D>;
%clear (FLOAT_T* array, int size);

%apply (FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* array, int size)};
%template(FLOAT_T_DESCRIPTOR_3D_LatticeUnSerializer) plb::NumPy2Lattice3D<FLOAT_T, plb::descriptors::DESCRIPTOR_3D>;
%clear (FLOAT_T* array, int size);

