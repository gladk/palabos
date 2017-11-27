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

public class MultiBlockLattice2D {
   private jlabos.double_D2Q9Descriptor_PlbMultiBlockLattice2D block;
   private Descriptor desc;
   private jlabos.Box2D domain;

   public MultiBlockLattice2D(jlabos.Box2D domain, Descriptor desc, BGK dynamic){
      this.desc = desc;
      this.domain = domain;
      if(desc == Descriptor.D2Q9){
         block = jlabos.double_d2q9.double_D2Q9Descriptor_generateMultiBlockLattice(domain, dynamic.generateDynamic2D());
      }else{
         System.out.println("Not yet implemented");
      }
   }

   public void defineDynamics(jlabos.MultiNTensorField2DInt mask, jlabos.BounceBack dynamic){
      jlabos.double_d2q9.double_D2Q9Descriptor_m_plbDefineDynamics(block, mask.getObject(), domain, dynamic.generateDynamic2D());	
   }

   public jlabos.double_D2Q9Descriptor_PlbMultiBlockLattice2D getObject(){
      return this.block;
   }

   public void setBoundaryVelocity(double[] velocity, jlabos.Box2D domain){
      jlabos.double_d2q9.double_D2Q9Descriptor_plbSetBoundaryVelocity(block, velocity, domain);
   }

   public void setBoundaryVelocity(jlabos.MultiNTensorField2DDouble mask, jlabos.Box2D domain){
      jlabos.double_d2q9.double_D2Q9Descriptor_plbSetBoundaryVelocity(block, mask.getObject(), domain);
   }

   public void initializeAtEquilibrium(jlabos.MultiNTensorField2DDouble pressure, jlabos.MultiNTensorField2DDouble velocity, jlabos.Box2D domain){
      jlabos.double_d2q9.double_D2Q9Descriptor_plbInitializeAtEquilibrium(block, pressure.getObject(), velocity.getObject(), domain);
   }

   public void initializeAtEquilibrium(double size, double[] val, jlabos.Box2D domain){
      jlabos.double_d2q9.double_D2Q9Descriptor_plbInitializeAtEquilibrium(block, size, val, domain);
   }

   public jlabos.MultiNTensorField2DDouble velocityNorm(jlabos.Box2D domain){
      return new MultiNTensorField2DDouble(jlabos.double_d2q9.double_D2Q9Descriptor_velocityNorm(block, domain));
   }

   public double averageEnergy(jlabos.Box2D domain){
      return jlabos.double_d2q9.double_D2Q9Descriptor_averageEnergy(block, domain);
   }

   public double averageEnergy(){
      return jlabos.double_d2q9.double_D2Q9Descriptor_averageEnergy(block, domain);
   }

   public jlabos.MultiNTensorField2DInt meshGrid(int comp){
      jlabos.PlbMultiNTensorField2DInt x = jlabos.double_d2q9.double_D2Q9Descriptor_generateIntNTensorField2D(block, domain, 1);
      jlabos.int_block.int_plbSetToCoordinate(x, domain, comp);
      return new jlabos.MultiNTensorField2DInt(x);
   }

   public jlabos.double_D2Q9Descriptor_PlbOuterBoxBC getOuterBoxBC(){
      return  jlabos.double_d2q9.double_D2Q9Descriptor_regularizedBC();
   }

   public void fromArray(double[] data, jlabos.Box2D domain){
      jlabos.double_D2Q9Descriptor_LatticeUnSerializer unserializer = new jlabos.double_D2Q9Descriptor_LatticeUnSerializer(block, domain);
      unserializer.execute(data, data.length);
   }
   

   public double[] toArray(jlabos.Box2D domain){
      jlabos.double_D2Q9Descriptor_LatticeSerializer serializer = new double_D2Q9Descriptor_LatticeSerializer(block, domain);
      double[] myArray = new double[serializer.getSize()];
      serializer.execute(myArray, serializer.getSize());
      return myArray;
   }

   public void collideAndStream(){
      if(desc == Descriptor.D2Q9){
         block.collideAndStream();
      }
   }
}
