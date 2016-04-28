/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef SMARTUQ_TAYLOR_FUNCTIONS_H
#define SMARTUQ_TAYLOR_FUNCTIONS_H

#include <iostream>
#include <limits>
#include "taylor.h"

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
taylor_polynomial<T> sin(const taylor_polynomial<T> &other);

template <class T>
/**
 * @brief cos overloaded cos function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function cos in a polynomial
 */
taylor_polynomial<T> cos(const taylor_polynomial<T> &other);

template <class T>
/**
 * @brief tan overloaded tan function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function tan in a polynomial
 */
taylor_polynomial<T> tan(const taylor_polynomial<T> &other);

template <class T>
/**
 * @brief atan overloaded sin function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function atan in a polynomial
 */
taylor_polynomial<T> atan(const taylor_polynomial<T> &other);


// OTHERS
template <class T>
/**
 * @brief exp overloaded exp function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function exp in a polynomial
 */
taylor_polynomial<T> exp(const taylor_polynomial<T> &other);

template <class T>
/**
 * @brief tan overloaded tanh function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function tanh in a polynomial
 */
taylor_polynomial<T> tanh(const taylor_polynomial<T> &other);

template <class T>
/**
 * @brief sqrt overloaded sqrt function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function sqrt in a polynomial
 */
taylor_polynomial<T> sqrt(const taylor_polynomial<T> &other);

template <class T>
/**
 * @brief log overloaded log function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @return the evaluation of the function log in a polynomial
 */
taylor_polynomial<T> log(const taylor_polynomial<T> &other);

template <class T>
/**
 * @brief pow overloaded pow function (evaluated in a polynomial value)
 * @param other polynomial for evaluation
 * @param exponent exponent value
 * @return the evaluation of the function pow in a polynomial
 */
taylor_polynomial<T> pow(const taylor_polynomial<T> &other, const int &exponent);



#endif // SMARTUQ_TAYLOR_FUNCTIONS_H
