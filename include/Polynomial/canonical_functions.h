#ifndef CANONICAL_FUNCTIONS_H
#define CANONICAL_FUNCTIONS_H


#include <iostream>
#include <limits>
#include "canonical.h"

using namespace std;
using namespace smart;
using namespace polynomial;

//TRIGONOMETRIC FUNCTIONS
template <class T>
canonical_polynomial<T> sin(const canonical_polynomial<T> &other);

template <class T>
canonical_polynomial<T> cos(const canonical_polynomial<T> &other);

template <class T>
canonical_polynomial<T> tan(const canonical_polynomial<T> &other);

template <class T>
canonical_polynomial<T> cot(const canonical_polynomial<T> &other);

template <class T>
canonical_polynomial<T> asin(const canonical_polynomial<T> &other);

template <class T>
canonical_polynomial<T> acos(const canonical_polynomial<T> &other);

template <class T>
canonical_polynomial<T> atan(const canonical_polynomial<T> &other);

template <class T>
canonical_polynomial<T> acot(const canonical_polynomial<T> &other);

// OTHERS
template <class T>
canonical_polynomial<T> exp(const canonical_polynomial<T> &other);

template <class T>
canonical_polynomial<T> sqrt(const canonical_polynomial<T> &other);

template <class T>
canonical_polynomial<T> log(const canonical_polynomial<T> &other);

template <class T>
canonical_polynomial<T> log10(const canonical_polynomial<T> &other);

template <class T>
canonical_polynomial<T> pow(const canonical_polynomial<T> &other, const int &exponent);

template <class T>
canonical_polynomial<T> pow(const canonical_polynomial<T> &other, const double &exponent);



#endif // CANONICAL_FUNCTIONS_H
