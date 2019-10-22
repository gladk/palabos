/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * Main author: Yann Sagon <yann.sagon@unige.ch>
 */

package jlabos;

public class MultiBlockLattice3D {
   private jlabos.double_D3Q19Descriptor_PlbMultiBlockLattice3D block;
   private Descriptor desc;
   private jlabos.Box3D domain;

   public MultiBlockLattice3D(jlabos.Box3D domain, Descriptor desc, BGK dynamic){
      this.desc = desc;
      this.domain = domain;
      if(desc == Descriptor.D3Q19){	
         block = jlabos.double_d3q19.double_D3Q19Descriptor_generateMultiBlockLattice(domain, dynamic.generateDynamic3D());
      }else{
         System.out.println("Not yet implemented");
      }
   }	

   public void defineDynamics(jlabos.MultiNTensorField3DInt mask, jlabos.BounceBack dynamic){
      jlabos.double_d3q19.double_D3Q19Descriptor_m_plbDefineDynamics(block, mask.getObject(), domain, dynamic.generateDynamic3D());
   }

   public jlabos.double_D3Q19Descriptor_PlbMultiBlockLattice3D getObject(){
      return this.block;
   }

   public jlabos.Box3D getDomain(){
      return this.domain;
   }

   public jlabos.MultiNTensorField3DDouble getPopulations(jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_d3q19.double_D3Q19Descriptor_populations(block, domain));
   }

   public void setBoundaryVelocity(double[] velocity, jlabos.Box3D domain){
      jlabos.double_d3q19.double_D3Q19Descriptor_plbSetBoundaryVelocity(block, velocity, domain);
   }

   public  void setBoundaryDensity(double rho, jlabos.Box3D domain){
      jlabos.double_d3q19.double_D3Q19Descriptor_plbSetBoundaryDensity(block, domain, rho); 
   }

   public void initializeAtEquilibrium(double size, double[] val, jlabos.Box3D domain){
      jlabos.double_d3q19.double_D3Q19Descriptor_plbInitializeAtEquilibrium(block, size, val, domain);
   }

   public jlabos.MultiNTensorField3DDouble velocityNorm(jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_d3q19.double_D3Q19Descriptor_velocityNorm(block, domain));
   }

   public jlabos.MultiNTensorField3DDouble strainRateFromStress(){
      return strainRateFromStress(domain);
   }

   public jlabos.MultiNTensorField3DDouble strainRateFromStress(jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_d3q19.double_D3Q19Descriptor_strainRateFromStress(block, domain));
   }

   public double averageEnergy(jlabos.Box3D domain){
      return jlabos.double_d3q19.double_D3Q19Descriptor_averageEnergy(block, domain);
   }

   public jlabos.double_D3Q19Descriptor_PlbOuterBoxBC getOuterBoxBC(){
      return  jlabos.double_d3q19.double_D3Q19Descriptor_regularizedBC();
   }

   public void fromArray(double[] data, jlabos.Box3D domain){
      jlabos.double_D3Q19Descriptor_LatticeUnSerializer unserializer = new jlabos.double_D3Q19Descriptor_LatticeUnSerializer(block, domain);
      unserializer.execute(data, data.length);
   } 

   public jlabos.MultiNTensorField3DDouble velocity(){
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_d3q19.double_D3Q19Descriptor_velocity(block, domain));
   }


   public jlabos.MultiNTensorField3DDouble velocityComponent (jlabos.MultiNTensorField3DInt mask, int comp) {
      return velocityComponent (mask, mask.getDomain(), comp);
   }

   public jlabos.MultiNTensorField3DDouble velocityComponent (jlabos.MultiNTensorField3DInt mask, jlabos.Box3D domain, int comp) {
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_d3q19.double_D3Q19Descriptor_m_velocityComponent(block, mask.getObject(), domain, comp));
   }

   public double averageEnergy(){
      return jlabos.double_d3q19.double_D3Q19Descriptor_averageEnergy(block, domain);
   }

   public jlabos.MultiNTensorField3DInt meshGrid(int comp){
      jlabos.PlbMultiNTensorField3DInt x = jlabos.double_d3q19.double_D3Q19Descriptor_generateIntNTensorField3D(block, domain, 1);
      jlabos.int_block.int_plbSetToCoordinate(x, domain, comp);
      return new jlabos.MultiNTensorField3DInt(x);
   }

   public double[] toArray(jlabos.Box3D domain){
      jlabos.double_D3Q19Descriptor_LatticeSerializer serializer = new double_D3Q19Descriptor_LatticeSerializer(block, domain);
      double[] myArray = new double[serializer.getSize()];
      serializer.execute(myArray, serializer.getSize());
      return myArray;
   }

   public void collideAndStream(){
      if(desc == Descriptor.D3Q19){
         block.collideAndStream();
      }
   }
}
