/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "../../include/Polynomial/taylor_functions.h"

using namespace smartuq;
using namespace polynomial;

/************************************************/
/*                  SIN                         */
/************************************************/
template <class T>
taylor_polynomial<T> sin(const taylor_polynomial<T> &other){

    return taylor_polynomial<T>(1,1);
}
template class taylor_polynomial<double>
sin(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
sin(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
sin(const taylor_polynomial<long double> &);


/************************************************/
/*                  COS                         */
/************************************************/
template <class T>
taylor_polynomial<T> cos(const taylor_polynomial<T> &other){

    return taylor_polynomial<T>(1,1);

}
template class taylor_polynomial<double>
cos(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
cos(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
cos(const taylor_polynomial<long double> &);

/************************************************/
/*                  TAN                         */
/************************************************/
template <class T>
taylor_polynomial<T> tan(const taylor_polynomial<T> &other){

    return taylor_polynomial<T>(1,1);

}
template class taylor_polynomial<double>
tan(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
tan(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
tan(const taylor_polynomial<long double> &);


/************************************************/
/*                  ASIN                        */
/************************************************/
template <class T>
taylor_polynomial<T> asin(const taylor_polynomial<T> &other){

    return taylor_polynomial<T>(1,1);

}
template class taylor_polynomial<double>
asin(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
asin(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
asin(const taylor_polynomial<long double> &);

/************************************************/
/*                  ACOS                        */
/************************************************/
template <class T>
taylor_polynomial<T> acos(const taylor_polynomial<T> &other){

    return taylor_polynomial<T>(1,1);

}
template class taylor_polynomial<double>
acos(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
acos(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
acos(const taylor_polynomial<long double> &);

/************************************************/
/*                  ATAN                        */
/************************************************/
template <class T>
taylor_polynomial<T> atan(const taylor_polynomial<T> &other){

    return taylor_polynomial<T>(1,1);

}
template class taylor_polynomial<double>
atan(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
atan(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
atan(const taylor_polynomial<long double> &);

// OTHERS
//EXPONENTIAL FUNCTION
template <class T>
/************************************************/
/*                  EXP                         */
/************************************************/
taylor_polynomial<T> exp(const taylor_polynomial<T> &other){
         int nvar =  other.get_nvar();
         int degree = other.get_degree();

         std::vector <T> coeffs = other.get_coeffs(); // p = c + n(x)
         T c = coeffs[0];

         taylor_polynomial<T> n(nvar,degree);
         coeffs[0] = 0.0;
         n.set_coeffs(coeffs);

         taylor_polynomial<T> res(nvar,degree, (T) exp(c));
         taylor_polynomial<T> term = res;

         for (int i=1; i<=degree; i++){
             term /= (T) i;
             term *= n;
             res+= term;
         }

         return res;
}
template class taylor_polynomial<double>
exp(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
exp(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
exp(const taylor_polynomial<long double> &);


/************************************************/
/*                  SQRT                        */
/************************************************/
template <class T>
taylor_polynomial<T> sqrt(const taylor_polynomial<T> &other){

         int nvar =  other.get_nvar();
         int degree = other.get_degree();

         std::vector <T> coeffs = other.get_coeffs(); // p = c + n(x)
         T c = coeffs[0];
         if (c<0) {
             std::cout<<"Error: negative argument in sqrt"<<std::endl;
             exit(EXIT_FAILURE);
         }

         taylor_polynomial<T> times(nvar,degree); // times = n/c
         coeffs[0] = 0.0;
         times.set_coeffs(coeffs);
         times/=c;

         taylor_polynomial<T> res(nvar,degree, (T) sqrt(c));
         taylor_polynomial<T> term = res;

         T exponent = 0.5;
         for (int i=1; i<=degree; i++){
             term *= exponent/i;
             term *= times;
             res+= term;
             exponent -= 1;
         }

         return res;

}
template class taylor_polynomial<double>
sqrt(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
sqrt(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
sqrt(const taylor_polynomial<long double> &);

/************************************************/
/*                  LOG                         */
/************************************************/
template <class T>
taylor_polynomial<T> log(const taylor_polynomial<T> &other){

    return taylor_polynomial<T>(1,1);

}
template class taylor_polynomial<double>
log(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
log(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
log(const taylor_polynomial<long double> &);

/************************************************/
/*                  LOG10                       */
/************************************************/
template <class T>
taylor_polynomial<T> log10(const taylor_polynomial<T> &other){

    return taylor_polynomial<T>(1,1);

}
template class taylor_polynomial<double>
log10(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
log10(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
log10(const taylor_polynomial<long double> &);

/************************************************/
/*                  POW                         */
/************************************************/
template <class T>
taylor_polynomial<T> pow(const taylor_polynomial<T> &other, const int &exponent){
    if(exponent<=1){
        smart_throw("Pow function with integer exponent, integer must be > 1");
    }
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    taylor_polynomial<T> res(nvar,degree);

    res = other;
    for(int i=1; i<exponent; i++)
        res *= other;

    return res;
}
template class taylor_polynomial<double>
pow(const taylor_polynomial<double> &, const int &);
template class taylor_polynomial<float>
pow(const taylor_polynomial<float> &, const int &);
template class taylor_polynomial<long double>
pow(const taylor_polynomial<long double> &, const int &);


