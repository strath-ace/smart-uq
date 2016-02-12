#ifndef CHEBYSHEV_FUNCTIONS_H
#define CHEBYSHEV_FUNCTIONS_H


#include <iostream>
#include <limits>
#include "chebyshev.h"

using namespace std;
using namespace smart;
using namespace polynomial;

//TRIGONOMETRIC FUNCTIONS
template <class T>
chebyshev_polynomial<T> sin(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> cos(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> tan(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> cot(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> asin(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> acos(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> atan(const chebyshev_polynomial<T> &other);

template <class T>
chebyshev_polynomial<T> acot(const chebyshev_polynomial<T> &other);

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

template <class T>
chebyshev_polynomial<T> pow(const chebyshev_polynomial<T> &other, const double &exponent);



#endif // CHEBYSHEV_FUNCTIONS_H
