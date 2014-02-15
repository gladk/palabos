namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * plbWrapper/lattice/latticeAnalysisWrapper3D.i
 */

%{
#include "JLABOS_ROOT/plbWrapper/lattice/latticeAnalysisWrapper3D.h"
#include "PALABOS_ROOT/src/dataProcessors/dataAnalysisWrapper3D.h"
%}

namespace plb {

/** Write verbatim the few functions which are taken from 
  * PALABOS_ROOT/dataProcessors/dataAnalysisWrapper3D.h
  */

template<typename T, class Descriptor> 
T computeAverageDensity(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain);

template<typename T, class Descriptor> 
T computeAverageEnergy(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain);


/** Write verbatim the few functions which are taken from 
  * JLABOS_ROOT/plbWrapper/lattice/latticeAnalysisWrapper3D.h
  */

%newobject pypal_computeDensity;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* pypal_computeDensity (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain );

%newobject maskedPypal_computeDensity;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* maskedPypal_computeDensity (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<int>& mask,
        Box3D domain );

%newobject pypal_computeKineticEnergy;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* pypal_computeKineticEnergy (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain );

%newobject maskedPypal_computeKineticEnergy;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* maskedPypal_computeKineticEnergy (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<int>& mask,
        Box3D domain );

%newobject pypal_computeVelocityNorm;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* pypal_computeVelocityNorm (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain );

%newobject maskedPypal_computeVelocityNorm;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* maskedPypal_computeVelocityNorm (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<int>& mask,
        Box3D domain );

%newobject pypal_computeVelocityComponent;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* pypal_computeVelocityComponent (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain, plint iComponent );

%newobject maskedPypal_computeVelocityComponent;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* maskedPypal_computeVelocityComponent (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<int>& mask,
        Box3D domain, plint iComponent );

%newobject pypal_computeVelocity;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* pypal_computeVelocity (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain );

%newobject maskedPypal_computeVelocity;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* maskedPypal_computeVelocity (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<int>& mask,
        Box3D domain );

%newobject pypal_computePiNeq;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* pypal_computePiNeq (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain );

%newobject maskedPypal_computePiNeq;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* maskedPypal_computePiNeq (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<int>& mask,
        Box3D domain );

%newobject pypal_computeShearStress;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* pypal_computeShearStress (
    MultiBlockLattice3D<T,Descriptor>& lattice,
    Box3D domain );

%newobject maskedPypal_computeShearStress;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* maskedPypal_computeShearStress (
    MultiBlockLattice3D<T,Descriptor>& lattice,
    MultiNTensorField3D<int>& mask,
    Box3D domain );

%newobject pypal_computeStrainRateFromStress;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* pypal_computeStrainRateFromStress (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain );

%newobject maskedPypal_computeStrainRateFromStress;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* maskedPypal_computeStrainRateFromStress (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<int>& mask,
        Box3D domain );

%newobject pypal_computePopulation;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* pypal_computePopulation (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain, plint iPop );

%newobject maskedPypal_computePopulation;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* maskedPypal_computePopulation (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<int>& mask,
        Box3D domain, plint iPop );

%newobject pypal_computePopulations;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* pypal_computePopulations (
        MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain );

%newobject maskedPypal_computePopulations;
template<typename T, class Descriptor> 
MultiNTensorField3D<T>* maskedPypal_computePopulations (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<int>& mask, Box3D domain );

%newobject compute_UPO_Rhs;
template<typename T, class Descriptor>
MultiNTensorField3D<T>* compute_UPO_Rhs (
                           MultiNTensorField3D<T>& lattice, Box3D domain, T omega );

%newobject masked_compute_UPO_Rhs;
template<typename T, class Descriptor>
MultiNTensorField3D<T>* masked_compute_UPO_Rhs (
                                 MultiNTensorField3D<T>& lattice,
                                 MultiNTensorField3D<int>& mask, Box3D domain, T omega );

%newobject UPO_ApplyJ;
template<typename T, class Descriptor>
MultiNTensorField3D<T>* UPO_ApplyJ (
                           MultiNTensorField3D<T>& f,
                           MultiNTensorField3D<T>& g,
                           Box3D domain, T omega );

%newobject masked_UPO_ApplyJ;
template<typename T, class Descriptor>
MultiNTensorField3D<T>* masked_UPO_ApplyJ (
                                 MultiNTensorField3D<T>& f,
                                 MultiNTensorField3D<T>& g,
                                 MultiNTensorField3D<int>& mask, Box3D domain, T omega );

%newobject compute_UPO_EnergyDerivative;
template<typename T, class Descriptor>
MultiNTensorField3D<T>* compute_UPO_EnergyDerivative (
                           MultiNTensorField3D<T>& lattice, Box3D domain );

%newobject masked_compute_UPO_EnergyDerivative;
template<typename T, class Descriptor>
MultiNTensorField3D<T>* masked_compute_UPO_EnergyDerivative (
                                 MultiNTensorField3D<T>& lattice,
                                 MultiNTensorField3D<int>& mask, Box3D domain );
}  // namespace plb


%template(FLOAT_T_DESCRIPTOR_3D_averageDensity) plb::computeAverageDensity<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_averageEnergy) plb::computeAverageEnergy<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_density) plb::pypal_computeDensity<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_kineticEnergy) plb::pypal_computeKineticEnergy<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_velocityNorm) plb::pypal_computeVelocityNorm<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_velocityComponent) plb::pypal_computeVelocityComponent<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_velocity) plb::pypal_computeVelocity<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_deviatoricStress) plb::pypal_computePiNeq<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_deviatoricStress) plb::pypal_computeShearStress<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_strainRateFromStress) plb::pypal_computeStrainRateFromStress<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_population) plb::pypal_computePopulation<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_populations) plb::pypal_computePopulations<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_UPO_Rhs) plb::compute_UPO_Rhs<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_UPO_ApplyJ) plb::UPO_ApplyJ<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_UPO_EnergyDerivative) plb::compute_UPO_EnergyDerivative<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;

%template(FLOAT_T_DESCRIPTOR_3D_m_density) plb::maskedPypal_computeDensity<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_kineticEnergy) plb::maskedPypal_computeKineticEnergy<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_velocityNorm) plb::maskedPypal_computeVelocityNorm<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_velocityComponent) plb::maskedPypal_computeVelocityComponent<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_velocity) plb::maskedPypal_computeVelocity<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_deviatoricStress) plb::maskedPypal_computePiNeq<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_deviatoricStress) plb::maskedPypal_computeShearStress<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_strainRateFromStress) plb::maskedPypal_computeStrainRateFromStress<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_population) plb::maskedPypal_computePopulation<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_populations) plb::maskedPypal_computePopulations<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_UPO_Rhs) plb::masked_compute_UPO_Rhs<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_UPO_ApplyJ) plb::masked_UPO_ApplyJ<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;
%template(FLOAT_T_DESCRIPTOR_3D_m_UPO_EnergyDerivative) plb::masked_compute_UPO_EnergyDerivative<FLOAT_T,plb::descriptors::DESCRIPTOR_3D>;

