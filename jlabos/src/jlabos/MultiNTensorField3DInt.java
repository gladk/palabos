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
import java.io.*;

public class MultiNTensorField3DInt {
   private jlabos.PlbMultiNTensorField3DInt block;
   private Descriptor desc;
   private jlabos.Box3D domain;

   public MultiNTensorField3DInt(jlabos.Box3D domain, int dim){
      System.gc();
      this.domain = domain;
      block = jlabos.int_block.int_generateMultiNTensorField3D(domain, dim);
   }

   public MultiNTensorField3DInt generateIntNTensorField3D(jlabos.Box3D intersection, int nDim){
      return new MultiNTensorField3DInt(jlabos.int_block.int_generateIntNTensorField3D(block, intersection, nDim));
   }

   public MultiNTensorField3DInt(jlabos.PlbMultiNTensorField3DInt block){
      this.block = block;
      domain = block.getBoundingBox();
   }

   public jlabos.MultiNTensorField3DDouble toDouble(){
      return new jlabos.MultiNTensorField3DDouble(jlabos.int_block.int_plbCopyConvertDouble(block, domain));
   }

   public jlabos.MultiNTensorField3DInt add(int[] val){
      return add(val, domain);
   }

   public jlabos.MultiNTensorField3DInt add(int[] val, jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DInt(jlabos.int_block.int_add(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField3DInt add(jlabos.MultiNTensorField3DInt val){
      return add(val, domain);
   }

   public jlabos.MultiNTensorField3DInt add(jlabos.MultiNTensorField3DInt val, jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DInt(jlabos.int_block.int_add(block, val.getObject(), domain));
   }

   public jlabos.MultiNTensorField3DInt subtract(int val){
      return subtract(new int[] {val}, domain);
   }

   public jlabos.MultiNTensorField3DInt subtract(int[] val){
      return subtract(val, domain);
   }

   public jlabos.MultiNTensorField3DInt subtract(int[] val, jlabos.Box3D domain){
      System.gc();
      return new jlabos.MultiNTensorField3DInt(jlabos.int_block.int_subtract(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField3DInt subtract(jlabos.MultiNTensorField3DInt val){
      return subtract(val, domain);
   }

   public jlabos.MultiNTensorField3DInt subtract(jlabos.MultiNTensorField3DInt val, jlabos.Box3D domain){
      System.gc();
      return new jlabos.MultiNTensorField3DInt(jlabos.int_block.int_subtract(block, val.getObject(), domain));
   }

   public jlabos.MultiNTensorField3DInt multiply(int[] val){
      return multiply(val, domain);
   }

   public jlabos.MultiNTensorField3DInt multiply(int[] val , jlabos.Box3D domain){
      System.gc();
      return new jlabos.MultiNTensorField3DInt(jlabos.int_block.int_multiply(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField3DInt multiply(jlabos.MultiNTensorField3DInt val){
      return multiply(val, domain);
   }

   public jlabos.MultiNTensorField3DInt multiply(jlabos.MultiNTensorField3DInt val, jlabos.Box3D domain){
      System.gc();
      return new jlabos.MultiNTensorField3DInt(jlabos.int_block.int_multiply(block, val.getObject(), domain));
   }

   public jlabos.MultiNTensorField3DInt toThePower(int val){
      return toThePower(new int[] {val}, domain);
   }

   public jlabos.MultiNTensorField3DInt toThePower(int[] val){
      return toThePower(val, domain);
   }

   public jlabos.MultiNTensorField3DInt toThePower(int[] val , jlabos.Box3D domain){
      System.gc();
      return new jlabos.MultiNTensorField3DInt(jlabos.int_block.int_toThePower(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField3DInt equals(int val, jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DInt(jlabos.int_block.int_equals(block, new int[] {val}, 1, domain));
   }

   public jlabos.MultiNTensorField3DInt equals(int[] val, jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DInt(jlabos.int_block.int_equals(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField3DInt lessThan(int val){
      return lessThan(new int[] {val}, domain);
   }

   public jlabos.MultiNTensorField3DInt lessThan(int[] val){
      return lessThan(val, domain);
   }

   public jlabos.MultiNTensorField3DInt lessThan(int[] val , jlabos.Box3D domain){
      System.gc();
      return new jlabos.MultiNTensorField3DInt(jlabos.int_block.int_lessThan(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField3DInt greaterThan(int val){
      return greaterThan(new int[] {val}, domain);
   }

   public jlabos.MultiNTensorField3DInt greaterThan(int[] val){
      return greaterThan(val, domain);
   }

   public jlabos.MultiNTensorField3DInt greaterThan(int[] val, jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DInt(jlabos.int_block.int_greaterThan(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField3DInt negate(){
      return negate(domain);
   }                                      

   public jlabos.MultiNTensorField3DInt negate(jlabos.Box3D domain){                                      
      return new jlabos.MultiNTensorField3DInt(jlabos.int_block.int_negate(block, domain));
   }

   public jlabos.MultiNTensorField3DInt meshGrid(int comp){
      jlabos.PlbMultiNTensorField3DInt x = jlabos.int_block.int_generateIntNTensorField3D(block, domain, 1);
      jlabos.int_block.int_plbSetToCoordinate(x, domain, comp);
      return new MultiNTensorField3DInt(x);
   }

   public jlabos.PlbMultiNTensorField3DInt getObject(){
      return this.block;
   }

   public jlabos.Box3D getDomain(){
      return this.domain;
   }
}
