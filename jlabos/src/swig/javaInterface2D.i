namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * blockLattice2d/numPyInterface2D.i
 */

%{
#include "PALABOS_ROOT/src/latticeBoltzmann/nearestNeighborLattices2D.h"
#include "PALABOS_ROOT/src/latticeBoltzmann/nearestNeighborLattices2D.hh"
#include "JLABOS_ROOT/plbWrapper/lattice/numPyInterface2D.h"
%}


namespace plb {
%include "arrays_java.i";
%apply double[] {double *};
template<typename T, class Descriptor>
class Lattice2NumPy2D {
public:
    Lattice2NumPy2D(MultiBlockLattice2D<T,Descriptor>& lattice_);
    Lattice2NumPy2D(MultiBlockLattice2D<T,Descriptor>& lattice_, Box2D const& domain);
    void execute(T* array, int size);
    int getSize() const;
};

template<typename T, class Descriptor>
class NumPy2Lattice2D {
public:
    NumPy2Lattice2D(MultiBlockLattice2D<T,Descriptor>& lattice_);
    NumPy2Lattice2D(MultiBlockLattice2D<T,Descriptor>& lattice_, Box2D const& domain);
    void execute(T* array, int size);
    int getSize() const;
};

}  // namespace plb

/*%apply (FLOAT_T* INPLACE_ARRAY1, int DIM1) {(FLOAT_T* array, int size)};*/
%template(FLOAT_T_DESCRIPTOR_2D_LatticeSerializer) plb::Lattice2NumPy2D<FLOAT_T, plb::descriptors::DESCRIPTOR_2D>;
/*%clear (FLOAT_T* array, int size);*/

/* %apply (FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* array, int size)};*/
%template(FLOAT_T_DESCRIPTOR_2D_LatticeUnSerializer) plb::NumPy2Lattice2D<FLOAT_T, plb::descriptors::DESCRIPTOR_2D>;
/*%clear (FLOAT_T* array, int size);*/
