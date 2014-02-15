namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * blockLattice2d/multiBlockLattice2D.i
 */

%{
#include "PALABOS_ROOT/src/latticeBoltzmann/nearestNeighborLattices2D.h"
#include "PALABOS_ROOT/src/latticeBoltzmann/nearestNeighborLattices2D.hh"
#include "PALABOS_ROOT/src/multiBlock/multiBlockLattice2D.h"
#include "JLABOS_ROOT/plbWrapper/lattice/multiBlockInfo2D.h"
%}


namespace plb {

template<typename T, class Descriptor>
class MultiBlockLattice2D {
public:
    ~MultiBlockLattice2D();
    void getNx() const;
    void getNy() const;
    virtual void collide(Box2D domain);
    virtual void collide();
    virtual void stream(Box2D domain);
    virtual void stream();
    virtual void collideAndStream(Box2D domain);
    virtual void collideAndStream();
    virtual void incrementTime();
    void initialize();
private:
    MultiBlockLattice2D();
};

template<typename T, class Descriptor>
class MultiBlockLatticeInfo {
public:
    MultiBlockLatticeInfo(MultiBlockLattice2D<T,Descriptor> const& multiBlock);
    plint getNx() const;
    plint getNy() const;
    plint getNumBlocks() const;
    Box2D getSmallestBlock() const;
    Box2D getLargestBlock() const;
    plint getNumAllocatedCells() const;
};


}  // namespace plb

%template(FLOAT_T_DESCRIPTOR_2D_PlbMultiBlockLattice2D) plb::MultiBlockLattice2D<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_PlbMultiBlockLatticeInfo2D) plb::MultiBlockLatticeInfo<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;

