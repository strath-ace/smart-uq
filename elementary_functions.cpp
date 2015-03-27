#include "elementary_functions.h"


/************************************************/
/*                  SIN                         */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> sin(const Chebyshev_Polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Chebyshev_Polynomial<T> res(nvar,degree);

    //chebyshev expansion of sin in [a,b]
    std::vector<T> coeffs = other.get_coeffs();
    T range = 0.0;
    for(int i=0; i<coeffs.size(); i++)
        range += fabs(coeffs[i]);

    //approximate sin [-range,range]
    std::vector<T> cheb_sin = Chebyshev_Polynomial<T>::cheb_approximation(sin,-range,range);
    //univariate composition
    std::vector<Chebyshev_Polynomial<T> > base = Chebyshev_Polynomial<T>::evaluate_base(other,-range,range);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_sin[i];
    }

    return res;
}
template class Chebyshev_Polynomial<double>
sin(const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
sin(const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
sin(const Chebyshev_Polynomial<long double> &);


/************************************************/
/*                  COS                         */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> cos(const Chebyshev_Polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Chebyshev_Polynomial<T> res(nvar,degree);

    //chebyshev expansion of cos in [a,b]
    std::vector<T> coeffs = other.get_coeffs();
    T range = 0.0;
    for(int i=0; i<coeffs.size(); i++)
        range += fabs(coeffs[i]);

    //approximate cos [-range,range]
    std::vector<T> cheb_cos = Chebyshev_Polynomial<T>::cheb_approximation(cos,-range,range);
    //univariate composition
    std::vector<Chebyshev_Polynomial<T> > base = Chebyshev_Polynomial<T>::evaluate_base(other,-range,range);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_cos[i];
    }

    return res;
}
template class Chebyshev_Polynomial<double>
cos(const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
cos(const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
cos(const Chebyshev_Polynomial<long double> &);

/************************************************/
/*                  TAN                         */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> tan(const Chebyshev_Polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Chebyshev_Polynomial<T> res(nvar,degree);

    //chebyshev expansion of tan in [a,b]
    std::vector<T> coeffs = other.get_coeffs();
    T range = 0.0;
    for(int i=0; i<coeffs.size(); i++)
        range += fabs(coeffs[i]);

    //approximate tan [-range,range]
    std::vector<T> cheb_tan = Chebyshev_Polynomial<T>::cheb_approximation(tan,-range,range);
    //univariate composition
    std::vector<Chebyshev_Polynomial<T> > base = Chebyshev_Polynomial<T>::evaluate_base(other,-range,range);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_tan[i];
    }

    return res;
}
template class Chebyshev_Polynomial<double>
tan(const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
tan(const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
tan(const Chebyshev_Polynomial<long double> &);

/************************************************/
/*                  COT                         */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> cot(const Chebyshev_Polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;

}
template class Chebyshev_Polynomial<double>
cot(const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
cot(const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
cot(const Chebyshev_Polynomial<long double> &);

/************************************************/
/*                  ASIN                        */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> asin(const Chebyshev_Polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;

}
template class Chebyshev_Polynomial<double>
asin(const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
asin(const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
asin(const Chebyshev_Polynomial<long double> &);

/************************************************/
/*                  ACOS                        */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> acos(const Chebyshev_Polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;

}
template class Chebyshev_Polynomial<double>
acos(const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
acos(const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
acos(const Chebyshev_Polynomial<long double> &);

/************************************************/
/*                  ATAN                        */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> atan(const Chebyshev_Polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Chebyshev_Polynomial<T> res(nvar,degree);

    //chebyshev expansion of atan in [a,b]
    std::vector<T> coeffs = other.get_coeffs();
    T range = 0.0;
    for(int i=0; i<coeffs.size(); i++)
        range += fabs(coeffs[i]);

    //approximate atan [-range,range]
    std::vector<T> cheb_atan = Chebyshev_Polynomial<T>::cheb_approximation(atan,-range,range);
    //univariate composition
    std::vector<Chebyshev_Polynomial<T> > base = Chebyshev_Polynomial<T>::evaluate_base(other,-range,range);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_atan[i];
    }

    return res;
}
template class Chebyshev_Polynomial<double>
atan(const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
atan(const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
atan(const Chebyshev_Polynomial<long double> &);

/************************************************/
/*                  ACOT                        */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> acot(const Chebyshev_Polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;

}
template class Chebyshev_Polynomial<double>
acot(const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
acot(const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
acot(const Chebyshev_Polynomial<long double> &);



// OTHERS
//EXPONENTIAL FUNCTION
template <class T>
/************************************************/
/*                  EXP                         */
/************************************************/
Chebyshev_Polynomial<T> exp(const Chebyshev_Polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Chebyshev_Polynomial<T> res(nvar,degree);

    //chebyshev expansion of exp in [a,b]
    std::vector<T> coeffs = other.get_coeffs();
    T range = 0.0;
    for(int i=0; i<coeffs.size(); i++)
        range += fabs(coeffs[i]);

    //approximate exp [-range,range]
    std::vector<T> cheb_exp = Chebyshev_Polynomial<T>::cheb_approximation(exp,-range,range);
    //univariate composition
    std::vector<Chebyshev_Polynomial<T> > base = Chebyshev_Polynomial<T>::evaluate_base(other,-range,range);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_exp[i];
    }

    return res;
}
template class Chebyshev_Polynomial<double>
exp(const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
exp(const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
exp(const Chebyshev_Polynomial<long double> &);


/************************************************/
/*                  SQRT                        */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> sqrt(const Chebyshev_Polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;

}
template class Chebyshev_Polynomial<double>
sqrt(const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
sqrt(const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
sqrt(const Chebyshev_Polynomial<long double> &);

/************************************************/
/*                  LOG                         */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> log(const Chebyshev_Polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;

}
template class Chebyshev_Polynomial<double>
log(const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
log(const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
log(const Chebyshev_Polynomial<long double> &);

/************************************************/
/*                  LOG10                       */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> log10(const Chebyshev_Polynomial<T> &other){

}
template class Chebyshev_Polynomial<double>
log10(const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
log10(const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
log10(const Chebyshev_Polynomial<long double> &);

/************************************************/
/*                  POW                         */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> pow(const Chebyshev_Polynomial<T> &other, const int &exponent){
    if(exponent<=1){
        std::cout<<"pow function with integer exponent, integer must be > 1"<<std::endl;
        exit(EXIT_FAILURE);
    }
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Chebyshev_Polynomial<T> res(nvar,degree);

    res = other;
    for(int i=1; i<exponent; i++)
        res *= other;

    return res;
}
template class Chebyshev_Polynomial<double>
pow(const Chebyshev_Polynomial<double> &, const int &);
template class Chebyshev_Polynomial<float>
pow(const Chebyshev_Polynomial<float> &, const int &);
template class Chebyshev_Polynomial<long double>
pow(const Chebyshev_Polynomial<long double> &, const int &);



template <class T>
Chebyshev_Polynomial<T> pow(const Chebyshev_Polynomial<T> &other, const double &exponent){
std::cout<<"NOT IMPLEMENTED"<<std::endl;

}
template class Chebyshev_Polynomial<double>
pow(const Chebyshev_Polynomial<double> &, const double &);
template class Chebyshev_Polynomial<float>
pow(const Chebyshev_Polynomial<float> &, const double &);
template class Chebyshev_Polynomial<long double>
pow(const Chebyshev_Polynomial<long double> &, const double &);

