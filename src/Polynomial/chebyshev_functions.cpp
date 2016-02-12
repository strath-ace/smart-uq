#include "Polynomial/chebyshev_functions.h"

using namespace smart;
using namespace polynomial;

/************************************************/
/*                  SIN                         */
/************************************************/
template <class T>
chebyshev_polynomial<T> sin(const chebyshev_polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree);

    //chebyshev expansion of sin in [a,b]
    std::vector<T> range = other.get_range();

    //approximate sin [-range,range]
    std::vector<T> cheb_sin = chebyshev_polynomial<T>::cheb_approximation(sin,range[0],range[1]);
    //univariate composition
    std::vector<chebyshev_polynomial<T> > base = chebyshev_polynomial<T>::evaluate_base(other,range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_sin[i];
    }

    return res;
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
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree);

    //chebyshev expansion of cos in [a,b]
    std::vector<T> range = other.get_range();

    //approximate cos [-range,range]
    std::vector<T> cheb_cos = chebyshev_polynomial<T>::cheb_approximation(cos,range[0],range[1]);
    //univariate composition
    std::vector<chebyshev_polynomial<T> > base = chebyshev_polynomial<T>::evaluate_base(other,range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_cos[i];
    }

    return res;
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
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree);

    //chebyshev expansion of tan in [a,b]
    std::vector<T> range = other.get_range();

    //approximate tan [-range,range]
    std::vector<T> cheb_tan = chebyshev_polynomial<T>::cheb_approximation(tan,range[0],range[1]);
    //univariate composition
    std::vector<chebyshev_polynomial<T> > base = chebyshev_polynomial<T>::evaluate_base(other,range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_tan[i];
    }

    return res;
}
template class chebyshev_polynomial<double>
tan(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
tan(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
tan(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  COT                         */
/************************************************/
template <class T>
chebyshev_polynomial<T> cot(const chebyshev_polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;
    return chebyshev_polynomial<T>(0,0);

}
template class chebyshev_polynomial<double>
cot(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
cot(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
cot(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  ASIN                        */
/************************************************/
template <class T>
chebyshev_polynomial<T> asin(const chebyshev_polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;
    return chebyshev_polynomial<T>(0,0);
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
    std::cout<<"NOT IMPLEMENTED"<<std::endl;
    return chebyshev_polynomial<T>(0,0);
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
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree);

    //chebyshev expansion of atan in [a,b]
    std::vector<T> range = other.get_range();

    //approximate atan [-range,range]
    std::vector<T> cheb_atan = chebyshev_polynomial<T>::cheb_approximation(atan,range[0],range[1]);
    //univariate composition
    std::vector<chebyshev_polynomial<T> > base = chebyshev_polynomial<T>::evaluate_base(other,range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_atan[i];
    }

    return res;
}
template class chebyshev_polynomial<double>
atan(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
atan(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
atan(const chebyshev_polynomial<long double> &);

/************************************************/
/*                  ACOT                        */
/************************************************/
template <class T>
chebyshev_polynomial<T> acot(const chebyshev_polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;
    return chebyshev_polynomial<T>(0,0);
}
template class chebyshev_polynomial<double>
acot(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
acot(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
acot(const chebyshev_polynomial<long double> &);



// OTHERS
//EXPONENTIAL FUNCTION
template <class T>
/************************************************/
/*                  EXP                         */
/************************************************/
chebyshev_polynomial<T> exp(const chebyshev_polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree);

    //chebyshev expansion of exp in [a,b]
    std::vector<T> range = other.get_range();

    //approximate exp [-range,range]
    std::vector<T> cheb_exp = chebyshev_polynomial<T>::cheb_approximation(exp,range[0],range[1]);
    //univariate composition
    std::vector<chebyshev_polynomial<T> > base = chebyshev_polynomial<T>::evaluate_base(other,range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_exp[i];
    }

    return res;
}
template class chebyshev_polynomial<double>
exp(const chebyshev_polynomial<double> &);
template class chebyshev_polynomial<float>
exp(const chebyshev_polynomial<float> &);
template class chebyshev_polynomial<long double>
exp(const chebyshev_polynomial<long double> &);


/************************************************/
/*                  SQRT                        */
/************************************************/
template <class T>
chebyshev_polynomial<T> sqrt(const chebyshev_polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree);

    //chebyshev expansion of exp in [a,b]
    std::vector<T> range = other.get_range();

    //approximate exp [-range,range]
    std::vector<T> cheb_sqrt = chebyshev_polynomial<T>::cheb_approximation(sqrt,range[0],range[1]);
    //univariate composition
    std::vector<chebyshev_polynomial<T> > base = chebyshev_polynomial<T>::evaluate_base(other,range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_sqrt[i];
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
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree);

    //chebyshev expansion of exp in [a,b]
    std::vector<T> range = other.get_range();

    //approximate exp [-range,range]
    std::vector<T> cheb_log = chebyshev_polynomial<T>::cheb_approximation(log,range[0],range[1]);
    //univariate composition
    std::vector<chebyshev_polynomial<T> > base = chebyshev_polynomial<T>::evaluate_base(other,range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_log[i];
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
    std::cout<<"NOT IMPLEMENTED"<<std::endl;
        return chebyshev_polynomial<T>(0,0);
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
        std::cout<<"pow function with integer exponent, integer must be > 1"<<std::endl;
        exit(EXIT_FAILURE);
    }
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree);

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



template <class T>
chebyshev_polynomial<T> pow(const chebyshev_polynomial<T> &other, const double &exponent){
std::cout<<"NOT IMPLEMENTED"<<std::endl;
    return chebyshev_polynomial<T>(0,0);
}
template class chebyshev_polynomial<double>
pow(const chebyshev_polynomial<double> &, const double &);
template class chebyshev_polynomial<float>
pow(const chebyshev_polynomial<float> &, const double &);
template class chebyshev_polynomial<long double>
pow(const chebyshev_polynomial<long double> &, const double &);
