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

public class MultiNTensorField3DDouble {
   private jlabos.PlbMultiNTensorField3DDouble block;
   private Descriptor desc;
   private jlabos.Box3D domain;

   public MultiNTensorField3DDouble(jlabos.Box3D domain, int dim){
      System.gc();
      this.domain = domain;
      block = jlabos.double_block.double_generateMultiNTensorField3D(domain, dim);
   }

   public MultiNTensorField3DDouble(jlabos.PlbMultiNTensorField3DDouble block){
      this.block = block;
      domain = block.getBoundingBox();
   }

   public jlabos.MultiNTensorField3DDouble add(double[] val){
      return add(val, domain);
   }

   public jlabos.MultiNTensorField3DDouble add(double[] val, jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_block.double_add(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField3DDouble add(jlabos.MultiNTensorField3DDouble val){
      return add(val, domain);
   }

   public jlabos.MultiNTensorField3DDouble add(jlabos.MultiNTensorField3DDouble val, jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_block.double_add(block, val.getObject(), domain));
   }

   public jlabos.MultiNTensorField3DDouble subtract(double[] val){
      return subtract(val, domain);
   }

   public jlabos.MultiNTensorField3DDouble subtract(double[] val, jlabos.Box3D domain){
      System.gc();
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_block.double_subtract(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField3DDouble subtract(jlabos.MultiNTensorField3DDouble val){
      return subtract(val, domain);
   }

   public jlabos.MultiNTensorField3DDouble subtract(jlabos.MultiNTensorField3DDouble val, jlabos.Box3D domain){
      System.gc();
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_block.double_subtract(block, val.getObject(), domain));
   }

   public jlabos.MultiNTensorField3DDouble multiply(double[] val){
      return multiply(val, domain);
   }

   public jlabos.MultiNTensorField3DDouble multiply(double[] val , jlabos.Box3D domain){
      System.gc();
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_block.double_multiply(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField3DDouble multiply(jlabos.MultiNTensorField3DDouble val){
      return multiply(val, domain);
   }

   public jlabos.MultiNTensorField3DDouble multiply(jlabos.MultiNTensorField3DDouble val, jlabos.Box3D domain){
      System.gc();
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_block.double_multiply(block, val.getObject(), domain));
   }

   public jlabos.MultiNTensorField3DDouble toThePower(double[] val){
      return toThePower(val, domain);
   }

   public jlabos.MultiNTensorField3DDouble toThePower(double[] val , jlabos.Box3D domain){
      System.gc();
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_block.double_toThePower(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField3DInt lessThan(double[] val){
      return lessThan(val, domain);
   }

   public jlabos.MultiNTensorField3DInt lessThan(double[] val , jlabos.Box3D domain){
      System.gc();
      return new jlabos.MultiNTensorField3DInt(jlabos.double_block.double_lessThan(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField3DInt greaterThan(double[] val){
      return greaterThan(val, domain);
   }

   public jlabos.MultiNTensorField3DInt greaterThan(double[] val, jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DInt(jlabos.double_block.double_greaterThan(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField3DDouble strainRate(jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_block.double_plbComputeStrainRate(block, domain));
   }

   public jlabos.MultiNTensorField3DDouble symmetricTensorNorm(jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_block.double_plbComputeSymmetricTensorNorm(block, domain));
   }

   public void writeToFile(String filename){                                       
      jlabos.Mpi mpi = new jlabos.Mpi();                                       
      double[] myArray = toArray();                                          
      if(mpi.getRankMpi()==0){                                             
         try{                                                            
            // Create file                                             
            FileWriter fstream = new FileWriter(filename);           
            BufferedWriter out = new BufferedWriter(fstream);      
            for (int x=0; x<myArray.length; x++) {               
               out.write(Double.toString(myArray[x]));         
               out.newLine();                                
            }                                              
            out.close();                                 
         }catch (Exception e){//Catch exception if any 
            System.err.println("Error: " + e.getMessage());       
         }                                                      
      }                                                       
   }                                                        

   public jlabos.MultiNTensorField3DDouble computeStrainRate(){
      return computeStrainRate(domain);
   }

   public jlabos.MultiNTensorField3DDouble computeStrainRate(jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_block.double_plbComputeStrainRate(block, domain));
   }

   public jlabos.MultiNTensorField3DDouble computeSymmetricTensorNorm(jlabos.Box3D domain){
      return new jlabos.MultiNTensorField3DDouble(jlabos.double_block.double_plbComputeSymmetricTensorNorm(block, domain));
   }

   public double average(){
      double[] result = {0.};
      int size = 1;
      jlabos.double_block.double_plbAverage(block, domain, result, size);
      return result[0];
   }

   public double[][] to2DArray(){
      int nx = domain.getNx();
      int ny = domain.getNy();
      int nz = domain.getNz();
      //assert domain.getNy() == 1 : "Matrix must be [x][y][z] with y==1";
      double[] ser = toArray();

      int a = 0;
      int b = 0;

      if(nx == 1){
         a = ny;
         b = nz;
      }else if(ny == 1){
         a = nx;
         b = nz;
      }else if(nz == 1){
         a = nx;
         b = ny;
      }else{
         assert false : "bad dimensions";
      }

      double[][] myArray = new double[a][b];
      for(int x=0; x<a; x++){
         for(int z=0; z<b; z++){
            myArray[x][z] = ser[x*b+z];
         }
      }
      return myArray;
   }

   public double[] toArray(){                           
      jlabos.double_NTensorFieldSerializer3D mat = new jlabos.double_NTensorFieldSerializer3D(block);
      double[] myArray = new double[mat.getSize()];   
      mat.execute(myArray, mat.getSize());          
      return myArray;                              
   }

   public jlabos.PlbMultiNTensorField3DDouble getObject(){
      return this.block;
   }
}
