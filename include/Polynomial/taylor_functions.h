/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef TAYLOR_FUNCTIONS_H
#define TAYLOR_FUNCTIONS_H

#include <iostream>
#include <limits>
#include "taylor.h"

using namespace std;
using namespace smart;
using namespace polynomial;

//TRIGONOMETRIC FUNCTIONS
template <class T>
taylor_polynomial<T> sin(const taylor_polynomial<T> &other);

template <class T>
taylor_polynomial<T> cos(const taylor_polynomial<T> &other);

template <class T>
taylor_polynomial<T> tan(const taylor_polynomial<T> &other);

template <class T>
taylor_polynomial<T> asin(const taylor_polynomial<T> &other);

template <class T>
taylor_polynomial<T> acos(const taylor_polynomial<T> &other);

template <class T>
taylor_polynomial<T> atan(const taylor_polynomial<T> &other);

// OTHERS
template <class T>
taylor_polynomial<T> exp(const taylor_polynomial<T> &other);

template <class T>
taylor_polynomial<T> sqrt(const taylor_polynomial<T> &other);

template <class T>
taylor_polynomial<T> log(const taylor_polynomial<T> &other);

template <class T>
taylor_polynomial<T> log10(const taylor_polynomial<T> &other);

template <class T>
taylor_polynomial<T> pow(const taylor_polynomial<T> &other, const int &exponent);



#endif // TAYLOR_FUNCTIONS_H
