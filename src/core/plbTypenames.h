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

#ifndef PLB_TYPENAMES_H
#define PLB_TYPENAMES_H

#include <string>

namespace plb {

struct DynamicNativeType {
    virtual ~DynamicNativeType() { }
    virtual DynamicNativeType* clone() const =0;
    virtual int getTypeSize() const =0;
};

template <typename T>
struct NativeType : public DynamicNativeType {
    static char const* getName();
    virtual NativeType<T>* clone() const { return new NativeType<T>(*this); }
    virtual int getTypeSize() const { return sizeof(T); }
};

class NativeTypeConstructor {
public:
    NativeTypeConstructor(std::string typeName);
    ~NativeTypeConstructor();
    NativeTypeConstructor(NativeTypeConstructor const& rhs);
    NativeTypeConstructor& operator=(NativeTypeConstructor const& rhs);
    void swap(NativeTypeConstructor& rhs);
    int getTypeSize() const { return nativeType->getTypeSize(); }
private:
    DynamicNativeType* nativeType;
};

}  // namespace plb

#endif  // PLB_TYPENAMES_H
