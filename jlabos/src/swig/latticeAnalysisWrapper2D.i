namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * plbWrapper/lattice/latticeAnalysisWrapper2D.i
 */

%{
#include "JLABOS_ROOT/plbWrapper/lattice/latticeAnalysisWrapper2D.h"
#include "PALABOS_ROOT/src/dataProcessors/dataAnalysisWrapper2D.h"
%}

namespace plb {

/** Write verbatim the few functions which are taken from 
  * PALABOS_ROOT/dataProcessors/dataAnalysisWrapper2D.h
  */

template<typename T, class Descriptor> 
T computeAverageDensity(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain);

template<typename T, class Descriptor> 
T computeAverageEnergy(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain);


/** Write verbatim the few functions which are taken from 
  * JLABOS_ROOT/plbWrapper/lattice/latticeAnalysisWrapper2D.h
  */

%newobject pypal_computeDensity;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* pypal_computeDensity (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain );

%newobject maskedPypal_computeDensity;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* maskedPypal_computeDensity (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<int>& mask,
        Box2D domain );

%newobject pypal_computeKineticEnergy;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* pypal_computeKineticEnergy (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain );

%newobject maskedPypal_computeKineticEnergy;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* maskedPypal_computeKineticEnergy (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<int>& mask,
        Box2D domain );

%newobject pypal_computeVelocityNorm;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* pypal_computeVelocityNorm (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain );

%newobject maskedPypal_computeVelocityNorm;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* maskedPypal_computeVelocityNorm (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<int>& mask,
        Box2D domain );

%newobject pypal_computeVelocityComponent;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* pypal_computeVelocityComponent (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain, plint iComponent );

%newobject maskedPypal_computeVelocityComponent;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* maskedPypal_computeVelocityComponent (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<int>& mask,
        Box2D domain, plint iComponent );

%newobject pypal_computeVelocity;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* pypal_computeVelocity (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain );

%newobject maskedPypal_computeVelocity;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* maskedPypal_computeVelocity (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<int>& mask,
        Box2D domain );

%newobject pypal_computePiNeq;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* pypal_computePiNeq (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain );

%newobject maskedPypal_computePiNeq;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* maskedPypal_computePiNeq (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<int>& mask,
        Box2D domain );

%newobject pypal_computeShearStress;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* pypal_computeShearStress (
    MultiBlockLattice2D<T,Descriptor>& lattice,
    Box2D domain );

%newobject maskedPypal_computeShearStress;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* maskedPypal_computeShearStress (
    MultiBlockLattice2D<T,Descriptor>& lattice,
    MultiNTensorField2D<int>& mask,
    Box2D domain );

%newobject pypal_computeStrainRateFromStress;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* pypal_computeStrainRateFromStress (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain );

%newobject maskedPypal_computeStrainRateFromStress;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* maskedPypal_computeStrainRateFromStress (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<int>& mask,
        Box2D domain );

%newobject pypal_computePopulation;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* pypal_computePopulation (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain, plint iPop );

%newobject maskedPypal_computePopulation;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* maskedPypal_computePopulation (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<int>& mask,
        Box2D domain, plint iPop );

%newobject pypal_computePopulations;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* pypal_computePopulations (
        MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain );

%newobject maskedPypal_computePopulations;
template<typename T, class Descriptor> 
MultiNTensorField2D<T>* maskedPypal_computePopulations (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<int>& mask, Box2D domain );

}


%template(FLOAT_T_DESCRIPTOR_2D_averageDensity) plb::computeAverageDensity<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_averageEnergy) plb::computeAverageEnergy<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_density) plb::pypal_computeDensity<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_kineticEnergy) plb::pypal_computeKineticEnergy<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_velocityNorm) plb::pypal_computeVelocityNorm<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_velocityComponent) plb::pypal_computeVelocityComponent<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_velocity) plb::pypal_computeVelocity<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_deviatoricStress) plb::pypal_computePiNeq<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_deviatoricStress) plb::pypal_computeShearStress<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_strainRateFromStress) plb::pypal_computeStrainRateFromStress<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_population) plb::pypal_computePopulation<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_populations) plb::pypal_computePopulations<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;

%template(FLOAT_T_DESCRIPTOR_2D_m_density) plb::maskedPypal_computeDensity<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_m_kineticEnergy) plb::maskedPypal_computeKineticEnergy<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_m_velocityNorm) plb::maskedPypal_computeVelocityNorm<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_m_velocityComponent) plb::maskedPypal_computeVelocityComponent<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_m_velocity) plb::maskedPypal_computeVelocity<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_m_deviatoricStress) plb::maskedPypal_computePiNeq<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_m_deviatoricStress) plb::maskedPypal_computeShearStress<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_m_strainRateFromStress) plb::maskedPypal_computeStrainRateFromStress<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_m_population) plb::maskedPypal_computePopulation<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;
%template(FLOAT_T_DESCRIPTOR_2D_m_populations) plb::maskedPypal_computePopulations<FLOAT_T,plb::descriptors::DESCRIPTOR_2D>;

