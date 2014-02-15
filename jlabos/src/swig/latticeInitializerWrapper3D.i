namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * plbWrapper/lattice/latticeInitializerWrapper3D.i
 */

%{
#include "JLABOS_ROOT/plbWrapper/lattice/latticeInitializerWrapper3D.h"
#include "PALABOS_ROOT/src/dataProcessors/dataInitializerWrapper3D.h"
%}

namespace plb {

/* From
 * PALABOS_ROOT/dataProcessors/dataInitializerWrapper3D.h
 */

/*
As of now, Composite-Dynamics has not yet been wrapped.
template<typename T, class Descriptor>
void setCompositeDynamics( MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
                           CompositeDynamics<T,Descriptor>* compositeDynamics );
*/

template<typename T, class Descriptor>
void setBoundaryDensity(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain, T rho);

template<typename T, class Descriptor>
void setExternalScalar( MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
                        int whichScalar, T externalScalar );

/* From
 * JLABOS_ROOT/plbWrapper/lattice/latticeInitializerWrapper3D.h
 */

template<typename T, class Descriptor>
void pypalDefineDynamics(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
                         Dynamics<T,Descriptor>* dynamics);

template<typename T, class Descriptor>
void maskedDefineDynamics( MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
                           Box3D domain, Dynamics<T,Descriptor>* dynamics );

template<typename T, class Descriptor>
void setBoundaryVelocity( MultiBlockLattice3D<T,Descriptor>& lattice,
                          T* velocity, int numDimIs2, Box3D domain );

template<typename T, class Descriptor>
void setBoundaryVelocity( MultiBlockLattice3D<T,Descriptor>& lattice,
                          MultiNTensorField3D<T>& velocity,
                          Box3D domain );

template<typename T, class Descriptor>
void maskedSetBoundaryVelocity( MultiBlockLattice3D<T,Descriptor>& lattice,
                                MultiNTensorField3D<T>& velocity,
                                MultiNTensorField3D<int>& mask,
                                Box3D domain );

template<typename T, class Descriptor>
void initializeAtEquilibrium( MultiBlockLattice3D<T,Descriptor>& lattice,
                              T density, T* velocity, int numDimIs2, Box3D domain );

template<typename T, class Descriptor>
void initializeAtEquilibrium( MultiBlockLattice3D<T,Descriptor>& lattice,
                              MultiNTensorField3D<T>& density,
                              MultiNTensorField3D<T>& velocity,
                              Box3D domain );

template<typename T, class Descriptor>
void maskedInitializeAtEquilibrium( MultiBlockLattice3D<T,Descriptor>& lattice,
                                    MultiNTensorField3D<T>& density,
                                    MultiNTensorField3D<T>& velocity,
                                    MultiNTensorField3D<int>& mask,
                                    Box3D domain );

template<typename T, class Descriptor>
void setExternalVector( MultiBlockLattice3D<T,Descriptor>& lattice,
                        int vectorStartsAt, T* externalVector, int numDimIs2, Box3D domain );

template<typename T, class Descriptor>
void setPopulations( MultiBlockLattice3D<T,Descriptor>& lattice,
                     T* populations, int numDimIsQ, Box3D domain );

template<typename T, class Descriptor>
void setPopulations( MultiBlockLattice3D<T,Descriptor>& lattice,
                     MultiNTensorField3D<T>& populations,
                     Box3D domain );

template<typename T, class Descriptor>
void maskedSetPopulations( MultiBlockLattice3D<T,Descriptor>& lattice,
                           T* populations, int numDimIsQ,
                           MultiNTensorField3D<int>& mask, Box3D domain );

template<typename T, class Descriptor>
void maskedSetPopulations( MultiBlockLattice3D<T,Descriptor>& lattice,
                           MultiNTensorField3D<T>& populations,
                           MultiNTensorField3D<int>& mask,
                           Box3D domain );


}  // namespace plb

%include FLOAT_T_block_array.i

%apply(FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* externalVector, int numDimIs2)};
%apply(FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* velocity, int numDimIs2)};
%apply(FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* populations, int numDimIsQ)};

%template(FLOAT_T_DESCRIPTOR_3D_plbDefineDynamics) plb::pypalDefineDynamics<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_plbDefineDynamics) plb::maskedDefineDynamics<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_plbSetBoundaryVelocity) plb::setBoundaryVelocity<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_plbSetBoundaryVelocity) plb::maskedSetBoundaryVelocity<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_plbSetBoundaryDensity) plb::setBoundaryDensity<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_plbInitializeAtEquilibrium) plb::initializeAtEquilibrium<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_plbInitializeAtEquilibrium) plb::maskedInitializeAtEquilibrium<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
/*%template(plbSetCompositeDynamics) plb::setCompositeDynamics<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;*/
%template(FLOAT_T_DESCRIPTOR_3D_plbSetExternalScalar) plb::setExternalScalar<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_plbSetExternalVector) plb::setExternalVector<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_plbSetPopulations) plb::setPopulations<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_plbSetPopulations) plb::maskedSetPopulations<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
