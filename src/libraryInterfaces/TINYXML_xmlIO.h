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

/** \file
 * Input/Output in XML format -- header file.
 */

#ifndef XML_IO_H
#define XML_IO_H

#include <string>
#include <vector>
#include <map>
#include "tinyxml/tinyxml.h"
#include "core/array.h"

namespace plb {

class XMLreaderProxy;

class XMLreader {
public:
    XMLreader( std::string fName );
    XMLreader( std::vector<TiXmlNode*> pParentVect );
    ~XMLreader();
    void print(int indent) const;
    XMLreaderProxy operator[] (std::string name) const;
    XMLreaderProxy getElement(std::string name, plint id) const;
    std::string getName() const;
    std::string getText() const;
    std::string getText(plint id) const;
    plint getFirstId() const;
    std::string getFirstText() const;
    bool idExists(plint id) const;
    bool getNextId(plint& id) const;
    std::vector<XMLreader*> const& getChildren(plint id) const;
private:
    void mainProcessorIni(TiXmlNode* pParent);
    void mainProcessorIni(std::vector<TiXmlNode*> pParentVect);
    void slaveProcessorIni();
    XMLreader();
private:
    struct Data {
        std::string text;
        std::vector<XMLreader*> children;
    };
private: 
    std::string name;
    std::map<plint,Data> data_map;
    static XMLreader notFound;
};

class XMLreaderProxy {
public:
    XMLreaderProxy(XMLreader const* reader_);
    XMLreaderProxy(XMLreader const* reader_, plint id_);
    template <typename T> void read(T& values) const;
    template <typename T> bool readNoThrow(T& values) const;
    template <typename T> void read(std::vector<T>& values) const;
    template <typename T> bool readNoThrow(std::vector<T>& values) const;
    template <typename T, plint N> void read(Array<T,N>& values) const;
    template <typename T, plint N> bool readNoThrow(Array<T,N>& values) const;
    XMLreaderProxy operator[] (std::string name) const;
    XMLreaderProxy operator[] (plint newId) const;
    bool isValid() const;
    plint getId() const;
    XMLreaderProxy iterId() const;
    std::string getName() const;
    std::vector<XMLreader*> const& getChildren() const;
private:
    XMLreader const* reader;
    plint id;
};

class XMLwriter {
public:
    XMLwriter();
    ~XMLwriter();
    template<typename T> void set(T const& value, plint precision=-1);
    void setString(std::string const& value);
    template<typename T> void set(std::vector<T> const& values, plint precision=-1);
    template<typename T, int N> void set(Array<T,N> const& values, plint precision=-1);
    XMLwriter& operator[] (std::string name);
    XMLwriter& operator[] (plint id);
    template<typename ostrT>
    void toOutputStream(ostrT& ostr, int indent=0) const;
    void print(std::string fName) const;
private:
    XMLwriter( std::string name_ );
private:
    struct Data {
        std::string text;
        std::vector<XMLwriter*> children;
    };
private: 
    bool isDocument;
    std::string name;
    std::map<plint,Data> data_map;
    plint currentId;
};

}  // namespace plb

#endif  // XML_IO_H
