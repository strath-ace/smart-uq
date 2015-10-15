#ifndef ELEMENTARY_FUNCTIONS_H
#define ELEMENTARY_FUNCTIONS_H


#include <iostream>
#include <limits>
#include "chebyshev_polynomial.h"

using namespace std;
using namespace smart;
using namespace intrusive;

//DIRECT MULTIPLICATION
template <class T>
Chebyshev_Polynomial<T> direct_multiplication(const Chebyshev_Polynomial<T> &x0, const Chebyshev_Polynomial<T> &x1);

//TRIGONOMETRIC FUNCTIONS
template <class T>
Chebyshev_Polynomial<T> sin(const Chebyshev_Polynomial<T> &other);

template <class T>
Chebyshev_Polynomial<T> cos(const Chebyshev_Polynomial<T> &other);

template <class T>
Chebyshev_Polynomial<T> tan(const Chebyshev_Polynomial<T> &other);

template <class T>
Chebyshev_Polynomial<T> cot(const Chebyshev_Polynomial<T> &other);

template <class T>
Chebyshev_Polynomial<T> asin(const Chebyshev_Polynomial<T> &other);

template <class T>
Chebyshev_Polynomial<T> acos(const Chebyshev_Polynomial<T> &other);

template <class T>
Chebyshev_Polynomial<T> atan(const Chebyshev_Polynomial<T> &other);

template <class T>
Chebyshev_Polynomial<T> acot(const Chebyshev_Polynomial<T> &other);

// OTHERS
template <class T>
Chebyshev_Polynomial<T> exp(const Chebyshev_Polynomial<T> &other);

template <class T>
Chebyshev_Polynomial<T> sqrt(const Chebyshev_Polynomial<T> &other);

template <class T>
Chebyshev_Polynomial<T> log(const Chebyshev_Polynomial<T> &other);

template <class T>
Chebyshev_Polynomial<T> log10(const Chebyshev_Polynomial<T> &other);

template <class T>
Chebyshev_Polynomial<T> pow(const Chebyshev_Polynomial<T> &other, const int &exponent);

template <class T>
Chebyshev_Polynomial<T> pow(const Chebyshev_Polynomial<T> &other, const double &exponent);



#endif // ELEMENTARY_FUNCTIONS_H