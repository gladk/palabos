namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * plbWrapper/lattice/latticeInitializerWrapper2D.i
 */

%{
#include "JLABOS_ROOT/plbWrapper/lattice/latticeInitializerWrapper2D.h"
#include "PALABOS_ROOT/src/dataProcessors/dataInitializerWrapper2D.h"
%}

namespace plb {

/* From
 * PALABOS_ROOT/dataProcessors/dataInitializerWrapper2D.h
 */
/*
As of now, Composite-Dynamics has not yet been wrapped.
template<typename T, class Descriptor>
void setCompositeDynamics( MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain,
                           CompositeDynamics<T,Descriptor>* compositeDynamics );
*/

template<typename T, class Descriptor>
void setBoundaryDensity(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain, T rho);

template<typename T, class Descriptor>
void setExternalScalar( MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain,
                        int whichScalar, T externalScalar );

/* From
 * JLABOS_ROOT/plbWrapper/lattice/latticeInitializerWrapper2D.h
 */

template<typename T, class Descriptor>
void pypalDefineDynamics(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain,
                         Dynamics<T,Descriptor>* dynamics);

template<typename T, class Descriptor>
void maskedDefineDynamics( MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
                           Box2D domain, Dynamics<T,Descriptor>* dynamics );

template<typename T, class Descriptor>
void setBoundaryVelocity( MultiBlockLattice2D<T,Descriptor>& lattice,
                          T* velocity, int numDimIs2, Box2D domain );

template<typename T, class Descriptor>
void setBoundaryVelocity( MultiBlockLattice2D<T,Descriptor>& lattice,
                          MultiNTensorField2D<T>& velocity,
                          Box2D domain );

template<typename T, class Descriptor>
void maskedSetBoundaryVelocity( MultiBlockLattice2D<T,Descriptor>& lattice,
                                MultiNTensorField2D<T>& velocity,
                                MultiNTensorField2D<int>& mask,
                                Box2D domain );

template<typename T, class Descriptor>
void initializeAtEquilibrium( MultiBlockLattice2D<T,Descriptor>& lattice,
                              T density, T* velocity, int numDimIs2, Box2D domain );

template<typename T, class Descriptor>
void initializeAtEquilibrium( MultiBlockLattice2D<T,Descriptor>& lattice,
                              MultiNTensorField2D<T>& density,
                              MultiNTensorField2D<T>& velocity,
                              Box2D domain );

template<typename T, class Descriptor>
void maskedInitializeAtEquilibrium( MultiBlockLattice2D<T,Descriptor>& lattice,
                                    MultiNTensorField2D<T>& density,
                                    MultiNTensorField2D<T>& velocity,
                                    MultiNTensorField2D<int>& mask,
                                    Box2D domain );

template<typename T, class Descriptor>
void setExternalVector( MultiBlockLattice2D<T,Descriptor>& lattice,
                        int vectorStartsAt, T* externalVector, int numDimIs2, Box2D domain );

template<typename T, class Descriptor>
void setPopulations( MultiBlockLattice2D<T,Descriptor>& lattice,
                     T* populations, int numDimIsQ, Box2D domain );

template<typename T, class Descriptor>
void setPopulations( MultiBlockLattice2D<T,Descriptor>& lattice,
                     MultiNTensorField2D<T>& populations,
                     Box2D domain );

template<typename T, class Descriptor>
void maskedSetPopulations( MultiBlockLattice2D<T,Descriptor>& lattice,
                           T* populations, int numDimIsQ,
                           MultiNTensorField2D<int>& mask, Box2D domain );

template<typename T, class Descriptor>
void maskedSetPopulations( MultiBlockLattice2D<T,Descriptor>& lattice,
                           MultiNTensorField2D<T>& populations,
                           MultiNTensorField2D<int>& mask,
                           Box2D domain );


}  // namespace plb
/*
%include array.i
*/
%include "FLOAT_T_block_array.i";
%apply(FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* externalVector, int numDimIs2)};
%apply(FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* velocity, int numDimIs2)};
%apply(FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* populations, int numDimIsQ)};

%template(FLOAT_T_DESCRIPTOR_2D_plbDefineDynamics) plb::pypalDefineDynamics<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_m_plbDefineDynamics) plb::maskedDefineDynamics<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_plbSetBoundaryVelocity) plb::setBoundaryVelocity<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_m_plbSetBoundaryVelocity) plb::maskedSetBoundaryVelocity<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_plbSetBoundaryDensity) plb::setBoundaryDensity<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_plbInitializeAtEquilibrium) plb::initializeAtEquilibrium<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_m_plbInitializeAtEquilibrium) plb::maskedInitializeAtEquilibrium<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
/*%template(plbSetCompositeDynamics) plb::setCompositeDynamics<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;*/
%template(FLOAT_T_DESCRIPTOR_2D_plbSetExternalScalar) plb::setExternalScalar<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_plbSetExternalVector) plb::setExternalVector<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_plbSetPopulations) plb::setPopulations<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_m_plbSetPopulations) plb::maskedSetPopulations<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
