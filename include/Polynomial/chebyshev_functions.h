/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef SMARTUQ_CHEBYSHEV_FUNCTIONS_H
#define SMARTUQ_CHEBYSHEV_FUNCTIONS_H


#include <iostream>
#include <limits>
#include "chebyshev.h"

using namespace std;
using namespace smartuq;
using namespace polynomial;

//TRIGONOMETRIC FUNCTIONS
template <class T>
/**
 * @brief sin overloaded sin function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function sin in a polynomial
 */
chebyshev_polynomial<T> sin(const chebyshev_polynomial<T> &other);

template <class T>
/**
 * @brief cos overloaded cos function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function cos in a polynomial
 */
chebyshev_polynomial<T> cos(const chebyshev_polynomial<T> &other);

template <class T>
/**
 * @brief tan overloaded tan function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function tan in a polynomial
 */
chebyshev_polynomial<T> tan(const chebyshev_polynomial<T> &other);

template <class T>
/**
 * @brief asin overloaded asin function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function asin in a polynomial
 */
chebyshev_polynomial<T> asin(const chebyshev_polynomial<T> &other);

template <class T>
/**
 * @brief acos overloaded acos function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function acos in a polynomial
 */
chebyshev_polynomial<T> acos(const chebyshev_polynomial<T> &other);

template <class T>
/**
 * @brief atan overloaded atan function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function atan in a polynomial
 */
chebyshev_polynomial<T> atan(const chebyshev_polynomial<T> &other);


// OTHERS
template <class T>
/**
 * @brief exp overloaded exp function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function exp in a polynomial
 */
chebyshev_polynomial<T> exp(const chebyshev_polynomial<T> &other);

template <class T>
/**
 * @brief sqrt overloaded sqrt function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function sqrt in a polynomial
 */
chebyshev_polynomial<T> sqrt(const chebyshev_polynomial<T> &other);

template <class T>
/**
 * @brief log overloaded log function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function log in a polynomial
 */
chebyshev_polynomial<T> log(const chebyshev_polynomial<T> &other);

template <class T>
/**
 * @brief log10 overloaded log10 function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function log10 in a polynomial
 */
chebyshev_polynomial<T> log10(const chebyshev_polynomial<T> &other);

template <class T>
/**
 * @brief pow overloaded pow function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @param exponent exponent value
 * @return the evaluation of the function pow in a polynomial
 */
chebyshev_polynomial<T> pow(const chebyshev_polynomial<T> &other, const int &exponent);


#endif // SMARTUQ_CHEBYSHEV_FUNCTIONS_H
