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
/*                 ATAN2                        */
/************************************************/
template <class T>
taylor_polynomial<T> atan2(const taylor_polynomial<T> &y, const taylor_polynomial<T> &x){

    T tix = x.get_coeffs()[0];
    T tiy = y.get_coeffs()[0];

    T titheta= atan2(tiy,tix);

    T tixy = sqrt(tix*tix+tiy*tiy);
    T sinti = tiy/tixy;
    T costi = tix/tixy;

    taylor_polynomial<T> xx = costi*x+sinti*y;
    taylor_polynomial<T> yy = -sinti*x+costi*y;

    return titheta+atan(yy/xx);

}
template class taylor_polynomial<double>
atan2(const taylor_polynomial<double> &, const taylor_polynomial<double> &);
template class taylor_polynomial<float>
atan2(const taylor_polynomial<float> &, const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
atan2(const taylor_polynomial<long double> &, const taylor_polynomial<long double> &);

/************************************************/
/*                  SIN                         */
/************************************************/
template <class T>
taylor_polynomial<T> sin(const taylor_polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();

    std::vector <T> coeffs = other.get_coeffs(); // p = c + n(x)
    T c = coeffs[0];
    T f = 1;

    taylor_polynomial<T> n(nvar,degree);
    coeffs[0] = 0.0;
    n.set_coeffs(coeffs);

    taylor_polynomial<T> res(nvar,degree, (T) sin(c));
    taylor_polynomial<T> term=n;

    std::vector<taylor_polynomial<T> > powers;
    std::vector<T> factorial;

    res += term*cos(c);

    for (int i=2; i<=degree; i++){
        term *= n;
        f *= i;
        powers.push_back(term);
        factorial.push_back(f);
    }

    for (int i=1; i<=degree/2.0;i++){
        res += sin(c)*(pow(-1,i)/(T)factorial[2*i-2])*powers[2*i-2];
    }

    for (int i=1; i<=degree/2.0-1;i++){
        res += cos(c)*(pow(-1,i)/(T)factorial[2*i+1-2])*powers[2*i+1-2];
    }

    return res;
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
    int nvar =  other.get_nvar();
    int degree = other.get_degree();

    std::vector <T> coeffs = other.get_coeffs(); // p = c + n(x)
    T c = coeffs[0];
    T f = 1;

    taylor_polynomial<T> n(nvar,degree);
    coeffs[0] = 0.0;
    n.set_coeffs(coeffs);

    taylor_polynomial<T> res(nvar,degree, (T) cos(c));
    taylor_polynomial<T> term=n;

    std::vector<taylor_polynomial<T> > powers;
    std::vector<T> factorial;

    res -= term*sin(c);

    for (int i=2; i<=degree; i++){
        term *=n;
        f *= i;
        powers.push_back(term);
        factorial.push_back(f);
    }

    for (int i=1; i<=degree/2.0;i++){
        res += cos(c)*(pow(-1,i)/(T)factorial[2*i-2])*powers[2*i-2];
    }

    for (int i=1; i<=degree/2.0-1;i++){
        res -= sin(c)*(pow(-1,i)/(T)factorial[2*i+1-2])*powers[2*i+1-2];
    }

    return res;

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

    return sin(other)/cos(other); //is this proper?

}
template class taylor_polynomial<double>
tan(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
tan(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
tan(const taylor_polynomial<long double> &);

/************************************************/
/*                  ASIN                         */
/************************************************/
template <class T>
taylor_polynomial<T> asin(const taylor_polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();

    std::vector <T> coeffs = other.get_coeffs(); // p = c + n(x)
    T c = coeffs[0];

    taylor_polynomial<T> n(nvar,degree);
    coeffs[0] = 0.0;
    n.set_coeffs(coeffs);

    taylor_polynomial<T> res(nvar,degree, (T) asin(c));
    taylor_polynomial<T> term=n;

    std::vector<taylor_polynomial<T> > powers;
    std::vector<T> sequence;

    sequence.push_back(asin(c));
    sequence.push_back(1.0/sqrt(1.0-c*c));
    sequence.push_back(sequence[1]*c/(2.0*(1-c*c)));
    sequence.push_back(sequence[1]*(2.0*c*c+1.0)/(6.0*(1.0-c*c)*(1.0-c*c)));
    powers.push_back(term);
    term *= n;    
    powers.push_back(term);
    term *= n;    
    powers.push_back(term);

    for (int i=4; i<=degree; i++){
        term *= n;
        powers.push_back(term);
        sequence.push_back(((i-2)*(i-2)*sequence[i-2]/(i-1)+c*(2*i-3)*sequence[i-1])/(i*(1.0-c*c)));    
    }

    for (int i=1; i<=degree;i++){
        res += sequence[i]*powers[i-1];
    }


    return res;
}
template class taylor_polynomial<double>
asin(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
asin(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
asin(const taylor_polynomial<long double> &);

/************************************************/
/*                  ACOS                         */
/************************************************/
template <class T>
taylor_polynomial<T> acos(const taylor_polynomial<T> &other){

    return M_PI/2.0-asin(other);
}
template class taylor_polynomial<double>
acos(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
acos(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
acos(const taylor_polynomial<long double> &);

/************************************************/
/*                  ATAN                         */
/************************************************/
template <class T>
taylor_polynomial<T> atan(const taylor_polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();

    std::vector <T> coeffs = other.get_coeffs(); // p = c + n(x)
    T c = coeffs[0];

    taylor_polynomial<T> n(nvar,degree);
    coeffs[0] = 0.0;
    n.set_coeffs(coeffs);

    taylor_polynomial<T> res(nvar,degree, (T) atan(c));
    taylor_polynomial<T> term=n;

    std::vector<taylor_polynomial<T> > powers;
    std::vector<T> sequence;

    sequence.push_back(atan(c));
    sequence.push_back(1.0/(1.0+c*c));
    sequence.push_back(-sequence[1]*c/(1.0+c*c));
    sequence.push_back(sequence[1]*(3.0*c*c-1.0)/(3.0*(1.0+c*c)*(1.0+c*c)));
    powers.push_back(term);
    term *= n;    
    powers.push_back(term);
    term *= n;    
    powers.push_back(term);

    for (int i=4; i<=degree; i++){
        term *= n;
        powers.push_back(term);
        sequence.push_back(-((i-2)*sequence[i-2]+2.0*c*(i-1)*sequence[i-1])/(i*(1.0+c*c)));    
    }

    for (int i=1; i<=degree;i++){
        res += sequence[i]*powers[i-1];
    }


    return res;
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
/*                  TANH                         */
/************************************************/
template <class T>
taylor_polynomial<T> tanh(const taylor_polynomial<T> &other){

    return (exp(other)-exp(-other))/(exp(other)+exp(-other));

}
template class taylor_polynomial<double>
tanh(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
tanh(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
tanh(const taylor_polynomial<long double> &);


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

    int nvar =  other.get_nvar();
    int degree = other.get_degree();

    std::vector <T> coeffs = other.get_coeffs(); // p = c + n(x)
    T c = coeffs[0];


    taylor_polynomial<T> n(nvar,degree);
    coeffs[0] = 0.0;
    n.set_coeffs(coeffs);

    taylor_polynomial<T> res(nvar,degree, (T) log(c));
    taylor_polynomial<T> tmp(nvar,degree);
    tmp = n;

    if (c==0) {
        std::cout<<"Error: division by zero in log"<<std::endl;
        exit(EXIT_FAILURE);
    }

    res += tmp/c;
    for (int i=2; i<=degree; i++){
        tmp *= n;
        res += pow(-1,i+1)*tmp/((T) i*pow(c,i));
    }

    return res;

}
template class taylor_polynomial<double>
log(const taylor_polynomial<double> &);
template class taylor_polynomial<float>
log(const taylor_polynomial<float> &);
template class taylor_polynomial<long double>
log(const taylor_polynomial<long double> &);


/************************************************/
/*                  POW                         */
/************************************************/
template <class T>
taylor_polynomial<T> pow(const taylor_polynomial<T> &other, const double &exponent){

    // NOTE: right now when doing pow(p(x),-k) it approximates x^(-k) and composes with p. Maybe it should approximate x^k and compose with 1/p?

    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    taylor_polynomial<T> res(nvar,degree);

    // trivial case
    if (exponent == 0) res = 1.0;

    // natural exponent
    else if ( exponent == floor(fabs(exponent)) ){//fancier way of checking if integer without doing exact comparison of two doubles?
        res = other;
        for(int i=1; i<exponent; i++)
            res *= other;
    }
    
    //general case
    else{

        std::vector <T> coeffs = other.get_coeffs(); // p = c + n(x)
        T c = coeffs[0];
        taylor_polynomial<T> n(nvar,degree);
        coeffs[0] = 0.0;
        n.set_coeffs(coeffs);

        res = pow(c,exponent); //for now, it is not my problem if C<0 and exponent is not integer and your polynomial turns into nan
        T k = exponent;
        taylor_polynomial<T> term = res;

        for (int i=1; i<=degree; i++){
        term /= (T) i;
        term /= c;
        term *= k;
        term *= n;
        res  += term;
        k--;
        }
    }

    return res;
}
template class taylor_polynomial<double>
pow(const taylor_polynomial<double> &, const double &);
template class taylor_polynomial<float>
pow(const taylor_polynomial<float> &, const double &);
template class taylor_polynomial<long double>
pow(const taylor_polynomial<long double> &, const double &);


