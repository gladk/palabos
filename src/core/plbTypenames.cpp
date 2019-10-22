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

#include "parallelism/mpiManager.h"
#include "core/plbTypenames.h"
#include "core/plbComplex.h"
#include "core/util.h"
#include "core/runTimeDiagnostics.h"
#include <algorithm>

namespace plb {

template<>
char const* NativeType<char>::getName() {
    static const char name[] = "char";
    return name;
}

template<>
char const* NativeType<bool>::getName() {
    static const char name[] = "bool";
    return name;
}

template<>
char const* NativeType<int>::getName() {
    static const char name[] = "int";
    return name;
}

template<>
char const* NativeType<long>::getName() {
    static const char name[] = "long";
    return name;
}

template<>
char const* NativeType<long long>::getName() {
    static const char name[] = "long-long";
    return name;
}

template<>
char const* NativeType<float>::getName() {
    static const char name[] = "float";
    return name;
}

template<>
char const* NativeType<double>::getName() {
    static const char name[] = "double";
    return name;
}

template<>
char const* NativeType<long double>::getName() {
    static const char name[] = "long-double";
    return name;
}

#ifdef PLB_USE_GCC
template<>
char const* NativeType<__float128>::getName() {
    static const char name[] = "float-128";
    return name;
}
#endif

template<>
char const* NativeType<Complex<float> >::getName() {
    static const char name[] = "complex-float";
    return name;
}

template<>
char const* NativeType<Complex<double> >::getName() {
    static const char name[] = "complex-double";
    return name;
}

template<>
char const* NativeType<Complex<long double> >::getName() {
    static const char name[] = "complex-long-double";
    return name;
}

#ifdef PLB_USE_GCC
template<>
char const* NativeType<Complex<__float128> >::getName() {
    static const char name[] = "complex-float-128";
    return name;
}
#endif

NativeTypeConstructor::NativeTypeConstructor(std::string typeName)
{
    std::string tn(util::tolower(typeName));
    if (tn=="char") {
        nativeType = new NativeType<char>;
    }
    else if (tn=="bool") {
        nativeType = new NativeType<bool>;
    }
    else if (tn=="int") {
        nativeType = new NativeType<int>;
    }
    else if (tn=="long") {
        nativeType = new NativeType<long>;
    }
    else if (tn=="long-long") {
        nativeType = new NativeType<long long>;
    }
    else if (tn=="float") {
        nativeType = new NativeType<float>;
    }
    else if (tn=="double") {
        nativeType = new NativeType<double>;
    }
    else if (tn=="long-double") {
        nativeType = new NativeType<long double>;
    }
#ifdef PLB_USE_GCC
    else if (tn=="float-128") {
        nativeType = new NativeType<__float128>;
    }
#endif
    else if (tn=="complex-float") {
        nativeType = new NativeType<Complex<float> >;
    }
    else if (tn=="complex-double") {
        nativeType = new NativeType<Complex<double> >;
    }
    else if (tn=="complex-long-double") {
        nativeType = new NativeType<Complex<long double> >;
    }
#ifdef PLB_USE_GCC
    else if (tn=="complex-float-128") {
        nativeType = new NativeType<Complex<__float128> >;
    }
#endif
    else {
        plbLogicError("Type "+typeName+" is not recognized.");
    }
}

NativeTypeConstructor::~NativeTypeConstructor() {
    delete nativeType;
}

NativeTypeConstructor::NativeTypeConstructor(NativeTypeConstructor const& rhs)
    : nativeType( rhs.nativeType->clone() )
{ }

NativeTypeConstructor& NativeTypeConstructor::operator=(NativeTypeConstructor const& rhs)
{
    NativeTypeConstructor(rhs).swap(*this);
    return *this;
}

void NativeTypeConstructor::swap(NativeTypeConstructor& rhs) {
    std::swap(nativeType, rhs.nativeType);
}

}  // namespace plb
