namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * blockLattice3d/multiBlockLattice3D.i
 */

%{
#include "PALABOS_ROOT/src/latticeBoltzmann/nearestNeighborLattices3D.h"
#include "PALABOS_ROOT/src/latticeBoltzmann/nearestNeighborLattices3D.hh"
#include "PALABOS_ROOT/src/multiBlock/multiBlockLattice3D.h"
#include "JLABOS_ROOT/plbWrapper/lattice/multiBlockInfo3D.h"
%}


namespace plb {

template<typename T, class Descriptor>
class MultiBlockLattice3D {
public:
    ~MultiBlockLattice3D();
    void getNx() const;
    void getNy() const;
    void getNz() const;
    virtual void collide(Box3D domain);
    virtual void collide();
    virtual void stream(Box3D domain);
    virtual void stream();
    virtual void collideAndStream(Box3D domain);
    virtual void collideAndStream();
    virtual void incrementTime();
    void initialize();
private:
    MultiBlockLattice3D();
};

template<typename T, class Descriptor>
class MultiBlockLatticeInfo {
public:
    MultiBlockLatticeInfo(MultiBlockLattice3D<T,Descriptor> const& multiBlock);
    plint getNx() const;
    plint getNy() const;
    plint getNz() const;
    plint getNumBlocks() const;
    Box3D getSmallestBlock() const;
    Box3D getLargestBlock() const;
    plint getNumAllocatedCells() const;
};


}  // namespace plb

%template(FLOAT_T_DESCRIPTOR_3D_PlbMultiBlockLattice3D) plb::MultiBlockLattice3D<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_PlbMultiBlockLatticeInfo3D) plb::MultiBlockLatticeInfo<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;

