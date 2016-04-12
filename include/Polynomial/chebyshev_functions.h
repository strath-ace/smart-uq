/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef CHEBYSHEV_FUNCTIONS_H
#define CHEBYSHEV_FUNCTIONS_H


#include <iostream>
#include <limits>
#include "chebyshev.h"

using namespace std;
using namespace smartuq;
using namespace polynomial;

//TRIGONOMETRIC FUNCTIONS
template <class T>
chebyshev_polynomial<T> sin(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> cos(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> tan(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> asin(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> acos(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> atan(const chebyshev_polynomial<T> &other);

// OTHERS
template <class T>
chebyshev_polynomial<T> exp(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> sqrt(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> log(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> log10(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> pow(const chebyshev_polynomial<T> &other, const int &exponent);


#endif // CHEBYSHEV_FUNCTIONS_H
