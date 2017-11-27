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

#ifndef PLB_COMPLEX_HH
#define PLB_COMPLEX_HH

#include "core/plbComplex.h"
#include <cmath>

namespace plb {

template<typename T>
T Complex<T>::pi = (T)4*std::atan((T)1);

template<typename T>
Complex<T>::Complex()
    : Re(), Imag()
{ }

template<typename T>
Complex<T>::Complex(T Re_)
    : Re(Re_), Imag(T())
{ }

template<typename T>
Complex<T>::Complex(T Re_, T Imag_)
    : Re(Re_), Imag(Imag_)
{ }

template<typename T>
template<typename U>
Complex<T>::operator U() const {
    return (U)Re;
}

template<typename T>
T Complex<T>::real() const {
    return Re;
}

template<typename T>
T Complex<T>::imaginary() const {
    return Imag;
}

template<typename T>
T Complex<T>::modulus() const {
    return std::sqrt(sqrModulus());
}

template<typename T>
T Complex<T>::sqrModulus() const {
    return Re*Re + Imag*Imag;
}

template<typename T>
Complex<T> Complex<T>::conjugate() const {
    return Complex<T>(Re, -Imag);
}

template<typename T>
T Complex<T>::argument() const {
    if (Re>T()) {
        return std::atan(Imag/Re);
    }
    else if (Re<T()) {
        return pi + std::atan(Imag/Re);
    }
    else {
        return pi/(T)2;
    }
}

template<typename T>
Complex<T> Complex<T>::intpow(int n) const {
    T r_pow_n = std::pow((T)modulus(), (T)n);
    T phi = argument();
    return Complex<T> (
             r_pow_n*(std::cos(n*phi)),
             r_pow_n*(std::sin(n*phi)) );
}

template<typename T>
Complex<T>& Complex<T>::operator+=(Complex<T> const& rhs) {
    Re += rhs.Re;
    Imag += rhs.Imag;
    return *this;
}

template<typename T>
template<typename U>
Complex<T>& Complex<T>::operator+=(U rhs) {
    Re += (T)rhs;
    return *this;
}

template<typename T>
Complex<T>& Complex<T>::operator-=(Complex<T> const& rhs) {
    Re -= rhs.Re;
    Imag -= rhs.Imag;
    return *this;
}

template<typename T>
template<typename U>
Complex<T>& Complex<T>::operator-=(U rhs) {
    Re -= (T)rhs;
    return *this;
}
template<typename T>
Complex<T> Complex<T>::operator-() const
{
    return Complex<T>(-Re, -Imag);
}

template<typename T>
Complex<T>& Complex<T>::operator*=(Complex<T> const& rhs) {
    T tmpRe = Re*rhs.Re - Imag*rhs.Imag;
    Imag = Imag*rhs.Re + Re*rhs.Imag;
    Re = tmpRe;
    return *this;
}

template<typename T>
template<typename U>
Complex<T>& Complex<T>::operator*=(U rhs) {
    Re *= (T)rhs;
    Imag *= (T)rhs;
    return *this;
}

template<typename T>
Complex<T>& Complex<T>::operator/=(Complex<T> const& rhs) {
    T rhsNormSqr = rhs.sqrModulus();
    T tmpRe = (Re*rhs.Re + Imag*rhs.Imag)/rhsNormSqr;
    Imag = (Imag*rhs.Re - Re*rhs.Imag)/rhsNormSqr;
    Re = tmpRe;
    return *this;
}


template<typename T>
template<typename U>
Complex<T>& Complex<T>::operator/=(U rhs) {
    Re /= (T)rhs;
    Imag /= (T)rhs;
    return *this;
}


template<typename T>
Complex<T> operator+(Complex<T> const& arg1, Complex<T> const& arg2) {
    return Complex<T>(arg1) += arg2;
}

template<typename T, typename U>
Complex<T> operator+(Complex<T> const& arg1, U arg2) {
    return Complex<T>(arg1) += (T)arg2;
}

template<typename T, typename U>
Complex<U> operator+(T arg1, Complex<U> const& arg2) {
    return Complex<U>((U)arg1+arg2.real(), (U)arg2.imaginary());
}



template<typename T>
Complex<T> operator-(Complex<T> const& arg1, Complex<T> const& arg2) {
    return Complex<T>(arg1) -= arg2;
}

template<typename T, typename U>
Complex<T> operator-(Complex<T> const& arg1, U arg2) {
    return Complex<T>(arg1) -= (T)arg2;
}

template<typename T, typename U>
Complex<U> operator-(T arg1, Complex<U> const& arg2) {
    return Complex<T>((U)arg1-arg2.real(), arg2.imaginary());
}



template<typename T>
Complex<T> operator*(Complex<T> const& arg1, Complex<T> const& arg2) {
    return Complex<T>(arg1) *= arg2;
}

template<typename T, typename U>
Complex<T> operator*(Complex<T> const& arg1, U arg2) {
    return Complex<T>(arg1) *= (T)arg2;
}

template<typename T, typename U>
Complex<U> operator*(T arg1, Complex<U> const& arg2) {
    return Complex<U>((U)arg1*arg2.real(), (U)arg1*arg2.imaginary());
}



template<typename T>
Complex<T> operator/(Complex<T> const& arg1, Complex<T> const& arg2) {
    return Complex<T>(arg1) /= arg2;
}

template<typename T, typename U>
Complex<T> operator/(Complex<T> const& arg1, U arg2) {
    return Complex<T>(arg1) /= (T)arg2;
}

template<typename T, typename U>
Complex<U> operator/(T arg1, Complex<U> const& arg2) {
    return Complex<U>((U)arg1/arg2.real(), (U)arg1/arg2.imaginary());
}

}  // namespace plb

#endif  // PLB_COMPLEX_HH
