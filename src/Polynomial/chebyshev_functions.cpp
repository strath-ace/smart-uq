/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "../../include/Polynomial/chebyshev_functions.h"

using namespace smartuq;
using namespace polynomial;

/************************************************/
/*                  SIN                         */
/************************************************/
template <class T>
chebyshev_polynomial<T> atan2(const chebyshev_polynomial<T> &y, const chebyshev_polynomial<T> &x){
    std::vector<T> x_range=x.get_range();
    std::vector<T> y_range=y.get_range();
    // T pi = 3.141592653589793;

    // if (x_range[0]>=0){ //1st, 4th or (1st and 4th) quadrants
    //     chebyshev_polynomial<T> rxy2 = x*x+y*y;
    //     return asin(y/sqrt(rxy2));
    // }
    // else if (y_range[0]>=0){ //2nd or (1st and 2nd) quadrants
    //     chebyshev_polynomial<T> rxy2 = x*x+y*y;
    //     return acos(x/sqrt(rxy2));
    // }
    // else if (y_range[1]<=0){ //3rd or (3rd and 4th) quadrants
    //     chebyshev_polynomial<T> rxy2 = x*x+y*y;
    //     return -acos(x/sqrt(rxy2));
    // }
    // else if (x_range[1]<=0){ //(2nd and 3rd) quadrants
    //     chebyshev_polynomial<T> rxy2 = x*x+y*y;
    //     return -asin(y/sqrt(rxy2))+pi;
    // }
    // else{ //all quadrants
    //     //should do it by bivariate interpolation and minimizing probability of discontinuity
        T tix = x.get_coeffs()[0];
        T tiy = y.get_coeffs()[0];

        T titheta= atan2(tiy,tix);

        T tixy = sqrt(tix*tix+tiy*tiy);
        T sinti = tiy/tixy;
        T costi = tix/tixy;

        chebyshev_polynomial<T> xx = costi*x-sinti*y;
        chebyshev_polynomial<T> yy = sinti*x+costi*y;

        // cout << titheta << endl;

        // chebyshev_polynomial<T> rxy2 = xx*xx+yy*yy;
        return titheta+atan(yy/xx);

        // cout << "Y, range = [" << y_range[0] <<"    ,    " << y_range[1] << " ]"<< endl;
        // //cout << y << endl;
        // cout << "X, range = [" << x_range[0] <<"    ,    " << x_range[1] << " ]"<< endl;
        // //cout << x << endl;
        // cout << "YY, range = [" << yy.get_range()[0] <<"    ,    " << yy.get_range()[1] << " ]"<< endl;
        // //cout << y << endl;
        // cout << "Xx, range = [" << xx.get_range()[0] <<"    ,    " << xx.get_range()[1] << " ]"<< endl;
        // //cout << x << endl;
        // cout << "Rxy2, range = [" << rxy2.get_range()[0] <<"    ,    " << rxy2.get_range()[1] << " ]"<< endl;
        // //smart_throw("atan2: current implementation does not allow angle in more than 2 quadrants");

    // }

}
template class chebyshev_polynomial<double>
atan2(const chebyshev_polynomial<double> &, const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
atan2(const chebyshev_polynomial<float> &, const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
atan2(const chebyshev_polynomial<long double> &, const chebyshev_polynomial<long double> &);

/************************************************/
/*                  SIN                         */
/************************************************/
template <class T>
chebyshev_polynomial<T> sin(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(sin,other);

}
template class chebyshev_polynomial<double>
sin(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
sin(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
sin(const chebyshev_polynomial<long double> &);


/************************************************/
/*                  COS                         */
/************************************************/
template <class T>
chebyshev_polynomial<T> cos(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(cos,other);

}
template class chebyshev_polynomial<double>
cos(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
cos(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
cos(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  TAN                         */
/************************************************/
template <class T>
chebyshev_polynomial<T> tan(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(tan,other);

}
template class chebyshev_polynomial<double>
tan(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
tan(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
tan(const chebyshev_polynomial<long double> &);


/************************************************/
/*                  ASIN                        */
/************************************************/
template <class T>
chebyshev_polynomial<T> asin(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(asin,other);

}
template class chebyshev_polynomial<double>
asin(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
asin(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
asin(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  ACOS                        */
/************************************************/
template <class T>
chebyshev_polynomial<T> acos(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(acos,other);

}
template class chebyshev_polynomial<double>
acos(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
acos(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
acos(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  ATAN                        */
/************************************************/
template <class T>
chebyshev_polynomial<T> atan(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(atan,other);

}
template class chebyshev_polynomial<double>
atan(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
atan(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
atan(const chebyshev_polynomial<long double> &);

// OTHERS
//EXPONENTIAL FUNCTION
template <class T>
/************************************************/
/*                  EXP                         */
/************************************************/
chebyshev_polynomial<T> exp(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(exp,other);

}
template class chebyshev_polynomial<double>
exp(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
exp(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
exp(const chebyshev_polynomial<long double> &);

//TANGENT HYPERBOLIC FUNCTION
template <class T>
/************************************************/
/*                  TANH                         */
/************************************************/
chebyshev_polynomial<T> tanh(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(tanh,other);

}
template class chebyshev_polynomial<double>
tanh(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
tanh(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
tanh(const chebyshev_polynomial<long double> &);


/************************************************/
/*                  SQRT                        */
/************************************************/
template <class T>
chebyshev_polynomial<T> sqrt(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(sqrt0,other);

}
template class chebyshev_polynomial<double>
sqrt(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
sqrt(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
sqrt(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  LOG                         */
/************************************************/
template <class T>
chebyshev_polynomial<T> log(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(log,other);

}
template class chebyshev_polynomial<double>
log(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
log(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
log(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  LOG10                       */
/************************************************/
template <class T>
chebyshev_polynomial<T> log10(const chebyshev_polynomial<T> &other){

    return chebyshev_polynomial<T>::approximation(log10,other);

}
template class chebyshev_polynomial<double>
log10(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
log10(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
log10(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  POW                         */
/************************************************/
template <class T>
chebyshev_polynomial<T> pow(const chebyshev_polynomial<T> &other, const int &exponent){
    if(exponent<=1){
        smart_throw("Pow function with integer exponent, integer must be > 1");
    }
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree,other.is_monomial_base());

    res = other;
    for(int i=1; i<exponent; i++)
        res *= other;

    return res;
}
template class chebyshev_polynomial<double>
pow(const chebyshev_polynomial<double> &, const int &);
template class chebyshev_polynomial<float>
pow(const chebyshev_polynomial<float> &, const int &);
template class chebyshev_polynomial<long double>
pow(const chebyshev_polynomial<long double> &, const int &);

