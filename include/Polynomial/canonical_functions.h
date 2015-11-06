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
Canonical_Polynomial<T> sin(const Canonical_Polynomial<T> &other);

template <class T>
Canonical_Polynomial<T> cos(const Canonical_Polynomial<T> &other);

template <class T>
Canonical_Polynomial<T> tan(const Canonical_Polynomial<T> &other);

template <class T>
Canonical_Polynomial<T> cot(const Canonical_Polynomial<T> &other);

template <class T>
Canonical_Polynomial<T> asin(const Canonical_Polynomial<T> &other);

template <class T>
Canonical_Polynomial<T> acos(const Canonical_Polynomial<T> &other);

template <class T>
Canonical_Polynomial<T> atan(const Canonical_Polynomial<T> &other);

template <class T>
Canonical_Polynomial<T> acot(const Canonical_Polynomial<T> &other);

// OTHERS
template <class T>
Canonical_Polynomial<T> exp(const Canonical_Polynomial<T> &other);

template <class T>
Canonical_Polynomial<T> sqrt(const Canonical_Polynomial<T> &other);

template <class T>
Canonical_Polynomial<T> log(const Canonical_Polynomial<T> &other);

template <class T>
Canonical_Polynomial<T> log10(const Canonical_Polynomial<T> &other);

template <class T>
Canonical_Polynomial<T> pow(const Canonical_Polynomial<T> &other, const int &exponent);

template <class T>
Canonical_Polynomial<T> pow(const Canonical_Polynomial<T> &other, const double &exponent);



#endif // CANONICAL_FUNCTIONS_H
