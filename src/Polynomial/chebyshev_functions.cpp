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
/*                 ATAN2                        */
/************************************************/
template <class T>
chebyshev_polynomial<T> atan2(const chebyshev_polynomial<T> &y, const chebyshev_polynomial<T> &x){

    T tix = x.get_coeffs()[0];
    T tiy = y.get_coeffs()[0];

    T titheta= atan2(tiy,tix);

    T tixy = sqrt(tix*tix+tiy*tiy);
    T sinti = tiy/tixy;
    T costi = tix/tixy;

    chebyshev_polynomial<T> xx = costi*x+sinti*y;
    chebyshev_polynomial<T> yy = -sinti*x+costi*y;

    return titheta+atan(yy/xx);

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

    //SPAGHETTI PATCH FOR RANGE CONSTRAINING
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree, other.is_monomial_base());
    std::vector<T> range = other.get_range();
    range[0]=max(range[0],(T) -1.0);
    range[1]=min(range[1],(T) 1.0);
    std::vector<T> approx(degree+1);


   if (other.is_monomial_base()){

       //approximate sin in [a,b] with increased degree (two-step truncation to enhance precision)
       int deg_max = chebyshev_polynomial<T>::MAX_DEGREE;
       int deg = std::min((int) (degree*1.5+1), deg_max);
       std::vector<T> cheb_approx = chebyshev_polynomial<T>::approximation(asin,range[0],range[1],deg);

       // Translation to canonical basis, taking into acount deg+1 terms from cheb_approx but building a monom_approx of degree+1 terms.
       // Hence rewriting code instead of calling to_monomial(), to avoid the computation of worthless terms of order > degree.

       chebyshev_polynomial<T> monom_approx(1,degree,(T) cheb_approx[0], true);
       chebyshev_polynomial<T> x(1,degree,(int) 0,-1.0,1.0, true);
       chebyshev_polynomial<T> cheb_base1(1,degree,(int) 0,-1.0,1.0, true);
       chebyshev_polynomial<T> cheb_base2(1,degree, (T) 1.0, true);
       chebyshev_polynomial<T> cheb_base(1,degree, true);

       monom_approx+= cheb_approx[1]*x;
       for (int i=2;i<=deg;i++){
           cheb_base=2.0*x*cheb_base1-cheb_base2;
           cheb_base2=cheb_base1;
           cheb_base1=cheb_base;
           monom_approx+=cheb_approx[i]*cheb_base;
       }

       approx = monom_approx.get_coeffs();
   }
   else{
       //approximate in [a,b]
       approx = chebyshev_polynomial<T>::approximation(asin,range[0],range[1],degree);
   }

    //univariate composition
    std::vector<chebyshev_polynomial<T> > base = chebyshev_polynomial<T>::evaluate_base1D(other,range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*approx[i];
    }

    return res;
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

    //SPAGHETTI PATCH FOR RANGE CONSTRAINING
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree, other.is_monomial_base());
    std::vector<T> range = other.get_range();
    range[0]=max(range[0],(T) -1.0);
    range[1]=min(range[1],(T) 1.0);
    std::vector<T> approx(degree+1);


   if (other.is_monomial_base()){

       //approximate sin in [a,b] with increased degree (two-step truncation to enhance precision)
       int deg_max = chebyshev_polynomial<T>::MAX_DEGREE;
       int deg = std::min((int) (degree*1.5+1), deg_max);
       std::vector<T> cheb_approx = chebyshev_polynomial<T>::approximation(acos,range[0],range[1],deg);

       // Translation to canonical basis, taking into acount deg+1 terms from cheb_approx but building a monom_approx of degree+1 terms.
       // Hence rewriting code instead of calling to_monomial(), to avoid the computation of worthless terms of order > degree.

       chebyshev_polynomial<T> monom_approx(1,degree,(T) cheb_approx[0], true);
       chebyshev_polynomial<T> x(1,degree,(int) 0,-1.0,1.0, true);
       chebyshev_polynomial<T> cheb_base1(1,degree,(int) 0,-1.0,1.0, true);
       chebyshev_polynomial<T> cheb_base2(1,degree, (T) 1.0, true);
       chebyshev_polynomial<T> cheb_base(1,degree, true);

       monom_approx+= cheb_approx[1]*x;
       for (int i=2;i<=deg;i++){
           cheb_base=2.0*x*cheb_base1-cheb_base2;
           cheb_base2=cheb_base1;
           cheb_base1=cheb_base;
           monom_approx+=cheb_approx[i]*cheb_base;
       }

       approx = monom_approx.get_coeffs();
   }
   else{
       //approximate in [a,b]
       approx = chebyshev_polynomial<T>::approximation(acos,range[0],range[1],degree);
   }

    //univariate composition
    std::vector<chebyshev_polynomial<T> > base = chebyshev_polynomial<T>::evaluate_base1D(other,range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*approx[i];
    }

    return res;
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

    //return chebyshev_polynomial<T>::approximation(sqrt,other);
    //SPAGHETTI PATCH FOR RANGE CONSTRAINING
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree, other.is_monomial_base());
    std::vector<T> range = other.get_range();
    range[0]=max(range[0],(T) (ZERO*ZERO));
    std::vector<T> approx(degree+1);


   if (other.is_monomial_base()){

       //approximate sin in [a,b] with increased degree (two-step truncation to enhance precision)
       int deg_max = chebyshev_polynomial<T>::MAX_DEGREE;
       int deg = std::min((int) (degree*1.5+1), deg_max);
       std::vector<T> cheb_approx = chebyshev_polynomial<T>::approximation(sqrt,range[0],range[1],deg);

       // Translation to canonical basis, taking into acount deg+1 terms from cheb_approx but building a monom_approx of degree+1 terms.
       // Hence rewriting code instead of calling to_monomial(), to avoid the computation of worthless terms of order > degree.

       chebyshev_polynomial<T> monom_approx(1,degree,(T) cheb_approx[0], true);
       chebyshev_polynomial<T> x(1,degree,(int) 0,-1.0,1.0, true);
       chebyshev_polynomial<T> cheb_base1(1,degree,(int) 0,-1.0,1.0, true);
       chebyshev_polynomial<T> cheb_base2(1,degree, (T) 1.0, true);
       chebyshev_polynomial<T> cheb_base(1,degree, true);

       monom_approx+= cheb_approx[1]*x;
       for (int i=2;i<=deg;i++){
           cheb_base=2.0*x*cheb_base1-cheb_base2;
           cheb_base2=cheb_base1;
           cheb_base1=cheb_base;
           monom_approx+=cheb_approx[i]*cheb_base;
       }

       approx = monom_approx.get_coeffs();
   }
   else{
       //approximate in [a,b]
       approx = chebyshev_polynomial<T>::approximation(sqrt,range[0],range[1],degree);
   }

    //univariate composition
    std::vector<chebyshev_polynomial<T> > base = chebyshev_polynomial<T>::evaluate_base1D(other,range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*approx[i];
    }

    return res;


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

    // return chebyshev_polynomial<T>::approximation(log,other);

    //SPAGHETTI PATCH FOR RANGE CONSTRAINING

    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree, other.is_monomial_base());
    std::vector<T> range = other.get_range();
    range[0]=max(range[0],(T) (ZERO*ZERO)); //where to set the min. a?? ZERO? exp(-1/ZERO)??
    std::vector<T> approx(degree+1);


    if (other.is_monomial_base()){

        //approximate sin in [a,b] with increased degree (two-step truncation to enhance precision)
        int deg_max = chebyshev_polynomial<T>::MAX_DEGREE;
        int deg = std::min((int) (degree*1.5+1), deg_max);
        std::vector<T> cheb_approx = chebyshev_polynomial<T>::approximation(log,range[0],range[1],deg);

        // Translation to canonical basis, taking into acount deg+1 terms from cheb_approx but building a monom_approx of degree+1 terms.
        // Hence rewriting code instead of calling to_monomial(), to avoid the computation of worthless terms of order > degree.

        chebyshev_polynomial<T> monom_approx(1,degree,(T) cheb_approx[0], true);
        chebyshev_polynomial<T> x(1,degree,(int) 0,-1.0,1.0, true);
        chebyshev_polynomial<T> cheb_base1(1,degree,(int) 0,-1.0,1.0, true);
        chebyshev_polynomial<T> cheb_base2(1,degree, (T) 1.0, true);
        chebyshev_polynomial<T> cheb_base(1,degree, true);

        monom_approx+= cheb_approx[1]*x;
        for (int i=2;i<=deg;i++){
            cheb_base=2.0*x*cheb_base1-cheb_base2;
            cheb_base2=cheb_base1;
            cheb_base1=cheb_base;
            monom_approx+=cheb_approx[i]*cheb_base;
        }

        approx = monom_approx.get_coeffs();
    }
    else{
        //approximate in [a,b]
        approx = chebyshev_polynomial<T>::approximation(log,range[0],range[1],degree);
    }

    //univariate composition
    std::vector<chebyshev_polynomial<T> > base = chebyshev_polynomial<T>::evaluate_base1D(other,range[0],range[1]);
    for (int i=0; i<=degree; i++){
     res += base[i]*approx[i];
    }

    return res;
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
    //THIS GUY NEEDS RANGE CONSTRAINING TOO
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
chebyshev_polynomial<T> pow(const chebyshev_polynomial<T> &other, const double &exponent){

    // NOTE: right now when doing pow(p(x),-k) it approximates x^(-k) and composes with p. Maybe it should approximate x^k and compose with 1/p?

    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree,other.is_monomial_base());

    // trivial case
    if (exponent == 0) res = 1.0;

    // natural exponent
    else if ( exponent == floor(fabs(exponent)) ){//fancier way of checking if integer without doing exact comparison of two doubles?
        res = other;
        for(int i=1; i<exponent; i++)
            res *= other;
    }

    //general case
    else{  //SPAGHETTI PATCH, DO THIS PROPERLY! 
    
        std::vector<T> range = other.get_range();
        if (exponent != floor (exponent)) range[0]=max(range[0],(T) 0.0); // constrain range if exponent is not integer
        // Note: this guy above is 0.0. That's good if exponent>1 but maybe should be pow(ZERO, 1/exponent) for exponent<1 to avioid approximating in a point of infinite derivative...? 
        std::vector<T> approx(degree+1);


        if (other.is_monomial_base()){

       //approximate sin in [a,b] with increased degree (two-step truncation to enhance precision)
        int deg_max = chebyshev_polynomial<T>::MAX_DEGREE;
        int deg = std::min((int) (degree*1.5+1), deg_max);

            int n = chebyshev_polynomial<T>::MAX_DEGREE;
            std::vector<T> cheb_approx(deg+1), d(n+1);
            T fac;
            T pi = 3.141592653589793;
            T t;
            T total;
            T y;

            for (int k = 0; k <= n; k++)
            {
                t = cos(pi*(k+0.5)/(n+1)); //zeros Ti
                y = ((1.0+t)*range[1] + (1.0-t)*range[0])/2.0; //mapped zeros
                d[k] = pow(y,exponent); //evaluate function
            }

            //Interpolation in chebyshev basis -> MAX_DEGREE chebyshev nodes but only a deg polynomial.
            fac = 2.0/(n+1);
            for (int j = 0; j <= deg; j++)
            {
                total = 0.0;
                for (int k = 0; k <= n; k++)
                {
                    total = total+d[k]*cos( (pi*j)*( (k+ 0.5)/(n+1) ) );
                }
                cheb_approx[j] = fac*total;
            }

            cheb_approx[0] = cheb_approx[0]/2.0;

            // Translation to canonical basis, taking into acount deg+1 terms from cheb_approx but building a monom_approx of degree+1 terms.
            // Hence rewriting code instead of calling to_monomial(), to avoid the computation of worthless terms of order > degree.

            chebyshev_polynomial<T> monom_approx(1,degree,(T) cheb_approx[0], true);
            chebyshev_polynomial<T> x(1,degree,(int) 0,-1.0,1.0, true);
            chebyshev_polynomial<T> cheb_base1(1,degree,(int) 0,-1.0,1.0, true);
            chebyshev_polynomial<T> cheb_base2(1,degree, (T) 1.0, true);
            chebyshev_polynomial<T> cheb_base(1,degree, true);

            monom_approx+= cheb_approx[1]*x;
            for (int i=2;i<=deg;i++){
               cheb_base=2.0*x*cheb_base1-cheb_base2;
               cheb_base2=cheb_base1;
               cheb_base1=cheb_base;
               monom_approx+=cheb_approx[i]*cheb_base;
            }

           approx = monom_approx.get_coeffs();
        }
        else{
            //approximate in [a,b]
                int n = chebyshev_polynomial<T>::MAX_DEGREE;
                std::vector<T> cheb_approx(degree+1), d(n+1);
                T fac;
                T pi = 3.141592653589793;
                T t;
                T total;
                T y;

                for (int k = 0; k <= n; k++)
                {
                    t = cos(pi*(k+0.5)/(n+1)); //zeros Ti
                    y = ((1.0+t)*range[1] + (1.0-t)*range[0])/2.0; //mapped zeros
                    d[k] = pow(y,exponent); //evaluate function
                }

                //Interpolation in chebyshev basis -> MAX_DEGREE chebyshev nodes but only a deg polynomial.
                fac = 2.0/(n+1);
                for (int j = 0; j <= degree; j++)
                {
                    total = 0.0;
                    for (int k = 0; k <= n; k++)
                    {
                        total = total+d[k]*cos( (pi*j)*( (k+ 0.5)/(n+1) ) );
                    }
                    cheb_approx[j] = fac*total;
                }

                cheb_approx[0] = cheb_approx[0]/2.0;
        }

        //univariate composition
        std::vector<chebyshev_polynomial<T> > base = chebyshev_polynomial<T>::evaluate_base1D(other,range[0],range[1]);
        for (int i=0; i<=degree; i++){
            res += base[i]*approx[i];
        }
    }

    return res;
}

template class chebyshev_polynomial<double>
pow(const chebyshev_polynomial<double> &, const double &);
template class chebyshev_polynomial<float>
pow(const chebyshev_polynomial<float> &, const double &);
template class chebyshev_polynomial<long double>
pow(const chebyshev_polynomial<long double> &, const double &);

