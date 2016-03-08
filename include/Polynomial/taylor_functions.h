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
