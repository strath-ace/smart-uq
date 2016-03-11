/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "Polynomial/chebyshev_functions.h"

using namespace smart;
using namespace polynomial;

/************************************************/
/*                  SIN                         */
/************************************************/
template <class T>
chebyshev_polynomial<T> sin(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(sin,other);

}
template class chebyshev_polynomial<double>
sin(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
sin(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
sin(const chebyshev_polynomial<long double> &);


/************************************************/
/*                  COS                         */
/************************************************/
template <class T>
chebyshev_polynomial<T> cos(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(cos,other);

}
template class chebyshev_polynomial<double>
cos(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
cos(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
cos(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  TAN                         */
/************************************************/
template <class T>
chebyshev_polynomial<T> tan(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(tan,other);

}
template class chebyshev_polynomial<double>
tan(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
tan(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
tan(const chebyshev_polynomial<long double> &);


/************************************************/
/*                  ASIN                        */
/************************************************/
template <class T>
chebyshev_polynomial<T> asin(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(asin,other);

}
template class chebyshev_polynomial<double>
asin(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
asin(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
asin(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  ACOS                        */
/************************************************/
template <class T>
chebyshev_polynomial<T> acos(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(acos,other);

}
template class chebyshev_polynomial<double>
acos(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
acos(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
acos(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  ATAN                        */
/************************************************/
template <class T>
chebyshev_polynomial<T> atan(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(atan,other);

}
template class chebyshev_polynomial<double>
atan(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
atan(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
atan(const chebyshev_polynomial<long double> &);

// OTHERS
//EXPONENTIAL FUNCTION
template <class T>
/************************************************/
/*                  EXP                         */
/************************************************/
chebyshev_polynomial<T> exp(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(exp,other);

}
template class chebyshev_polynomial<double>
exp(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
exp(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
exp(const chebyshev_polynomial<long double> &);


/************************************************/
/*                  SQRT                        */
/************************************************/
template <class T>
chebyshev_polynomial<T> sqrt(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(sqrt,other);

}
template class chebyshev_polynomial<double>
sqrt(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
sqrt(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
sqrt(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  LOG                         */
/************************************************/
template <class T>
chebyshev_polynomial<T> log(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(log,other);

}
template class chebyshev_polynomial<double>
log(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
log(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
log(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  LOG10                       */
/************************************************/
template <class T>
chebyshev_polynomial<T> log10(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(log10,other);

}
template class chebyshev_polynomial<double>
log10(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
log10(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
log10(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  POW                         */
/************************************************/
template <class T>
chebyshev_polynomial<T> pow(const chebyshev_polynomial<T> &other, const int &exponent){
    if(exponent<=1){
        smart_throw("Pow function with integer exponent, integer must be > 1");
    }
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree,other.is_monomial_base());

    res = other;
    for(int i=1; i<exponent; i++)
        res *= other;

    return res;
}
template class chebyshev_polynomial<double>
pow(const chebyshev_polynomial<double> &, const int &);
template class chebyshev_polynomial<float>
pow(const chebyshev_polynomial<float> &, const int &);
template class chebyshev_polynomial<long double>
pow(const chebyshev_polynomial<long double> &, const int &);

