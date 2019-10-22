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

public class MultiNTensorField2DDouble {
   private jlabos.PlbMultiNTensorField2DDouble block;
   private Descriptor desc;
   private jlabos.Box2D domain;

   public MultiNTensorField2DDouble(jlabos.Box2D domain, int dim){
      this.domain = domain;
      block = jlabos.double_block.double_generateMultiNTensorField2D(domain, dim);
   }

   public MultiNTensorField2DDouble(jlabos.PlbMultiNTensorField2DDouble block){
      this.block = block;
      domain = block.getBoundingBox();
   }

   public jlabos.MultiNTensorField2DDouble add(double[] val){
      return add(val, domain);
   }

   public jlabos.MultiNTensorField2DDouble add(double[] val, jlabos.Box2D domain){
      return new jlabos.MultiNTensorField2DDouble(jlabos.double_block.double_add(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField2DDouble add(jlabos.MultiNTensorField2DDouble val){
      return add(val, domain);
   }

   public jlabos.MultiNTensorField2DDouble add(jlabos.MultiNTensorField2DDouble val, jlabos.Box2D domain){
      return new jlabos.MultiNTensorField2DDouble(jlabos.double_block.double_add(block, val.getObject(), domain));
   }

   public jlabos.MultiNTensorField2DDouble subtract(double[] val){
      return subtract(val, domain);
   }   

   public jlabos.MultiNTensorField2DDouble subtract(double[] val, jlabos.Box2D domain){
      return new jlabos.MultiNTensorField2DDouble(jlabos.double_block.double_subtract(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField2DDouble subtract(jlabos.MultiNTensorField2DDouble val){
      return subtract(val, domain);
   }

   public jlabos.MultiNTensorField2DDouble subtract(jlabos.MultiNTensorField2DDouble val, jlabos.Box2D domain){
      return new jlabos.MultiNTensorField2DDouble(jlabos.double_block.double_subtract(block, val.getObject(), domain));
   }

   public jlabos.MultiNTensorField2DDouble multiply(double[] val){
      return multiply(val, domain);
   }

   public jlabos.MultiNTensorField2DDouble multiply(double[] val , jlabos.Box2D domain){
      return new jlabos.MultiNTensorField2DDouble(jlabos.double_block.double_multiply(block, val, val.length, domain));
   }

   public jlabos.MultiNTensorField2DDouble multiply(jlabos.MultiNTensorField2DDouble val){
      return multiply(val, domain);
   }

   public jlabos.MultiNTensorField2DDouble multiply(jlabos.MultiNTensorField2DDouble val, jlabos.Box2D domain){
      return new jlabos.MultiNTensorField2DDouble(jlabos.double_block.double_multiply(block, val.getObject(), domain));
   }

   public jlabos.MultiNTensorField2DDouble toThePower(double[] val){
      return toThePower(val, domain);
   }

   public jlabos.MultiNTensorField2DDouble toThePower(double[] val , jlabos.Box2D domain){
      return new jlabos.MultiNTensorField2DDouble(jlabos.double_block.double_toThePower(block, val, val.length, domain));	
   }

   public jlabos.MultiNTensorField2DInt lessThan(double[] val){
      return lessThan(val, domain);
   }

   public jlabos.MultiNTensorField2DInt lessThan(double[] val , jlabos.Box2D domain){
      return new jlabos.MultiNTensorField2DInt(jlabos.double_block.double_lessThan(block, val, val.length, domain));	
   }

   public void setToConstant(jlabos.Box2D domain, double[] val){
      jlabos.double_block.double_plbSetToConstant(block, domain, val, val.length);
   }
   public void setComponent(int component, jlabos.MultiNTensorField2DDouble pois, jlabos.Box2D domain){
      jlabos.double_block.double_plbAssignComponent(block, component, pois.getObject(), domain);
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

   public double[][] to2DArray(){
      double[] ser = toArray();

      double[][] myArray = new double[domain.getNx()][domain.getNy()];
      for(int x=0; x<domain.getNx(); x++){
         for(int y=0; y<domain.getNy(); y++){
            myArray[x][y] = ser[x*domain.getNy()+y];
         }
      }
      return myArray;
   }   

   public double[] toArray(){
      jlabos.double_NTensorFieldSerializer2D mat = new jlabos.double_NTensorFieldSerializer2D(block);
      double[] myArray = new double[mat.getSize()];
      mat.execute(myArray, mat.getSize());
      return myArray;
   }

   public void fromArray(double[] data){
      jlabos.double_NTensorFieldUnSerializer2D unserializer = new jlabos.double_NTensorFieldUnSerializer2D(block);
      unserializer.execute(data, data.length);
        // util.assert_ndarray_equals(serial_array, self._domain, self.ndim_elements())
      //unserializer = self.NTensorFieldUnSerializer()
      //unserializer.execute(serial_array.astype(self.dtype()).flatten())
   }

   public jlabos.PlbMultiNTensorField2DDouble getObject(){
      return this.block;
   }
}
