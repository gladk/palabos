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

package jlabos;

import org.yaml.snakeyaml.Yaml;
import java.io.*;
import java.util.*;
public class Param {
   private Yaml yaml;
   private InputStream input;
   private Map<String, Object> object;
   private Map<String, Object> list1;

   public Param(String file) throws java.io.FileNotFoundException {
      this.input = new FileInputStream(new File(file));
      this.yaml = new Yaml();
      this.object = (Map<String, Object>) yaml.load(input);
   }

   public Integer getIntValue(List<String> list) throws java.io.FileNotFoundException {
      Integer toto = null;
      return getValue(toto, list);
   }

   public Double getDblValue(List<String> list) throws java.io.FileNotFoundException {
      Double toto = null;
      return getValue(toto, list);
   }

   public String getStrValue(List<String> list) throws java.io.FileNotFoundException {
      String toto = null;
      return getValue(toto, list);
   }

   private <T> T getValue(T dest, List<String> list) throws java.io.FileNotFoundException {
      Map<String, Object> list1 = object;
      for(int i = 0; i< list.size()-1; i++){
         try{
            list1 = (Map<String, Object>) list1.get(list.get(i));
         }catch(NullPointerException e){
            jlabos.Utils.pout("Error reading param: " + list.get(i));
            java.lang.System.exit(1);
         }
      }
      return(T)(list1.get(list.get(list.size()-1)));
   }
}
