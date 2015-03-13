#include "elementary_functions.h"

template <class T>
Chebyshev_Polynomial<T> sin(const Chebyshev_Polynomial<T> &other){

}

template <class T>
Chebyshev_Polynomial<T> cos(const Chebyshev_Polynomial<T> &other){

}

template <class T>
Chebyshev_Polynomial<T> tan(const Chebyshev_Polynomial<T> &other){

}

template <class T>
Chebyshev_Polynomial<T> cot(const Chebyshev_Polynomial<T> &other){

}

template <class T>
Chebyshev_Polynomial<T> asin(const Chebyshev_Polynomial<T> &other){

}

template <class T>
Chebyshev_Polynomial<T> acos(const Chebyshev_Polynomial<T> &other){

}

template <class T>
Chebyshev_Polynomial<T> atan(const Chebyshev_Polynomial<T> &other){

}

template <class T>
Chebyshev_Polynomial<T> acot(const Chebyshev_Polynomial<T> &other){

}


// OTHERS
//EXPONENTIAL FUNCTION
template <class T>
Chebyshev_Polynomial<T> exp(const Chebyshev_Polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Chebyshev_Polynomial<T> res(nvar,degree);

    //chebyshev expansion of exp in [-1,1]
    Chebyshev_Polynomial<T> f_exp(nvar,degree);
    std::vector<T> cheb_exp(Chebyshev_Polynomial<T>::MAX_DEGREE+1);
    cheb_exp[0] = 1.266065877752008;
    cheb_exp[1] = 1.130318207984970;
    cheb_exp[2] = 0.271495339534077;
    cheb_exp[3] = 0.044336849848664;
    cheb_exp[4] = 0.005474240442094;
    cheb_exp[5] = 0.000542926311914;
    cheb_exp[6] = 0.000044977322954;
    cheb_exp[7] = 0.000003198436462;
    cheb_exp[8] = 0.000000199212481;
    cheb_exp[9] = 0.000000011036772;
    cheb_exp[10] = 0.000000000550590;
    cheb_exp[11] = 0.000000000024980;
    cheb_exp[12] = 0.000000000001039;
    cheb_exp[13] = 0.000000000000040;
    cheb_exp[14] = 0.000000000000001;
    std::vector<T> cheb_exp_coeff(f_exp.get_coeffs().size());
    //fill f_c
    //expansion of univariate elementary functions in [-1,1]
    cheb_exp_coeff[0] = cheb_exp[0];
    for(int i=0; i<degree; i++){
        int elements = factorial(nvar+i)/(factorial(nvar)*factorial(i));
        cheb_exp_coeff[elements] = cheb_exp[i+1];
    }
    f_exp.set_coeffs(cheb_exp_coeff);

    //computing exp(tau(other))
    res = f_exp.f_composition(other);

    return res;
}
template class Chebyshev_Polynomial<double>
exp(const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
exp(const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
exp(const Chebyshev_Polynomial<long double> &);


template <class T>
Chebyshev_Polynomial<T> sqrt(const Chebyshev_Polynomial<T> &other){

}

template <class T>
Chebyshev_Polynomial<T> log(const Chebyshev_Polynomial<T> &other){

}

template <class T>
Chebyshev_Polynomial<T> log10(const Chebyshev_Polynomial<T> &other){

}

template <class T>
Chebyshev_Polynomial<T> pow(const Chebyshev_Polynomial<T> &other, const int exponent){

}

template <class T>
Chebyshev_Polynomial<T> pow(const Chebyshev_Polynomial<T> &other, const double &exponent){

}

template <class T>
Chebyshev_Polynomial<T> inv(const Chebyshev_Polynomial<T> &other){

}
