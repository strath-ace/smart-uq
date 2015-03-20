#ifndef ELEMENTARY_FUNCTIONS_H
#define ELEMENTARY_FUNCTIONS_H


#include <iostream>
#include <limits>
#include "chebyshev_polynomial.h"

using namespace std;

//COMPOSITION OF CHEBYSHEV POLYNOMIALS
//evaluation of a univariate chebyshev polynomial in a chebyshev polynomial
template <class T>
Chebyshev_Polynomial<T> composition(const std::vector<T> &coeffs, const Chebyshev_Polynomial<T> &other, const T &range = 1.0);

//chebyshev approximation of univariate function over [a,b]
template <class T>
std::vector<T> cheb_approximation(T (*f)(T x), const T a, const T b);

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

template <class T>
Chebyshev_Polynomial<T> inv(const Chebyshev_Polynomial<T> &other);


#endif // ELEMENTARY_FUNCTIONS_H
