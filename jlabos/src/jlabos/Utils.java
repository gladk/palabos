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
import indiji.mlplot.MLPlot;
import indiji.mlplot.MLPlot.Style;
import indiji.mlplot.MLPlot.Symbol;
import java.io.*;

public class Utils {

   public static <T> void pout(T s){
      jlabos.JlabosBase base = jlabos.JlabosBase.getSingletonObject();
      if(base.getRankMpi()==0){
         System.out.println(s);
      }
   }
   
   public static void saveAsJpg(jlabos.MultiNTensorField2DDouble snap, String filename){
      save(snap, filename, "jpg", 0);
   }

   public static void saveAsJpg(jlabos.MultiNTensorField2DDouble snap, String filename, int size){
      save(snap, filename, "jpg", size);
   }

   public static void saveAsSvg(jlabos.MultiNTensorField2DDouble snap, String filename){
      save(snap, filename, "svg", 0);
   }

   public static void saveAsSvg(jlabos.MultiNTensorField2DDouble snap, String filename, int size){
      save(snap, filename, "svg", size);
   }

   public static void saveAsJpg(jlabos.MultiNTensorField3DDouble snap, String filename){
      save(snap, filename, "jpg", 0);
   }

   public static void saveAsJpg(jlabos.MultiNTensorField3DDouble snap, String filename, int size){
      save(snap, filename, "jpg", size);
   }

   public static void saveAsSvg(jlabos.MultiNTensorField3DDouble snap, String filename){
      save(snap, filename, "svg", 0);
   }

   public static void saveAsSvg(jlabos.MultiNTensorField3DDouble snap, String filename, int size){
      save(snap, filename, "svg", size);
   }

   private static void save(jlabos.MultiNTensorField2DDouble snap, String filename, String ext, int size){
        internalSave(snap.to2DArray(), filename, ext, size);
   }

   private static void save(jlabos.MultiNTensorField3DDouble snap, String filename, String ext, int size){
        internalSave(snap.to2DArray(), filename, ext, size);
   }

   private static void internalSave(double[][] data, String filename, String ext, int size){
      jlabos.JlabosBase base = jlabos.JlabosBase.getSingletonObject();
      if(base.getRankMpi()==0){
         MLPlot p=new MLPlot();

         // Draw data
         p.imagesc(data);
         p.setxTickLabelRotation(45);
         p.setxLabelDist(30);

         // Customize Plot
         p.setTitle("Image-Map, rotated xTickLabels, one of many ColorMaps,..");

         // Save Vector/Bitmap Graphics
         if(size!=0){
            p.save(new File(filename + "." + ext),size);
         }else{
            p.save(new File(filename + "." + ext));
         }
      }

   }

}
