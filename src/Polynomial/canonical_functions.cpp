#include "Polynomial/canonical_functions.h"

using namespace smart;
using namespace polynomial;

/************************************************/
/*                  SIN                         */
/************************************************/
template <class T>
canonical_polynomial<T> sin(const canonical_polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    canonical_polynomial<T> res(nvar,degree);

    //Canonical expansion of sin in [a,b]
    std::vector<T> range = other.get_range();

    //approximate sin [-range,range]
    std::vector<T> cheb_sin = canonical_polynomial<T>::approximation_1d(sin,range[0],range[1],other.get_degree());
    //univariate composition
    std::vector<canonical_polynomial<T> > base = other.evaluate_base(range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_sin[i];
    }

    return res;
}
template class canonical_polynomial<double>
sin(const canonical_polynomial<double> &);
template class canonical_polynomial<float>
sin(const canonical_polynomial<float> &);
template class canonical_polynomial<long double>
sin(const canonical_polynomial<long double> &);


/************************************************/
/*                  COS                         */
/************************************************/
template <class T>
canonical_polynomial<T> cos(const canonical_polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    canonical_polynomial<T> res(nvar,degree);

    //Canonical expansion of cos in [a,b]
    std::vector<T> range = other.get_range();

    //approximate cos [-range,range]
    std::vector<T> cheb_cos = canonical_polynomial<T>::approximation_1d(cos,range[0],range[1],other.get_degree());
    //univariate composition
    std::vector<canonical_polynomial<T> > base = other.evaluate_base(range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_cos[i];
    }

    return res;
}
template class canonical_polynomial<double>
cos(const canonical_polynomial<double> &);
template class canonical_polynomial<float>
cos(const canonical_polynomial<float> &);
template class canonical_polynomial<long double>
cos(const canonical_polynomial<long double> &);

/************************************************/
/*                  TAN                         */
/************************************************/
template <class T>
canonical_polynomial<T> tan(const canonical_polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    canonical_polynomial<T> res(nvar,degree);

    //Canonical expansion of tan in [a,b]
    std::vector<T> range = other.get_range();

    //approximate tan [-range,range]
    std::vector<T> cheb_tan = canonical_polynomial<T>::approximation_1d(tan,range[0],range[1],other.get_degree());
    //univariate composition
    std::vector<canonical_polynomial<T> > base = other.evaluate_base(range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_tan[i];
    }

    return res;
}
template class canonical_polynomial<double>
tan(const canonical_polynomial<double> &);
template class canonical_polynomial<float>
tan(const canonical_polynomial<float> &);
template class canonical_polynomial<long double>
tan(const canonical_polynomial<long double> &);

/************************************************/
/*                  COT                         */
/************************************************/
template <class T>
canonical_polynomial<T> cot(const canonical_polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;
    return canonical_polynomial<T>(0,0);
}
template class canonical_polynomial<double>
cot(const canonical_polynomial<double> &);
template class canonical_polynomial<float>
cot(const canonical_polynomial<float> &);
template class canonical_polynomial<long double>
cot(const canonical_polynomial<long double> &);

/************************************************/
/*                  ASIN                        */
/************************************************/
template <class T>
canonical_polynomial<T> asin(const canonical_polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;
    return canonical_polynomial<T>(0,0);
}
template class canonical_polynomial<double>
asin(const canonical_polynomial<double> &);
template class canonical_polynomial<float>
asin(const canonical_polynomial<float> &);
template class canonical_polynomial<long double>
asin(const canonical_polynomial<long double> &);

/************************************************/
/*                  ACOS                        */
/************************************************/
template <class T>
canonical_polynomial<T> acos(const canonical_polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;
    return canonical_polynomial<T>(0,0);
}
template class canonical_polynomial<double>
acos(const canonical_polynomial<double> &);
template class canonical_polynomial<float>
acos(const canonical_polynomial<float> &);
template class canonical_polynomial<long double>
acos(const canonical_polynomial<long double> &);

/************************************************/
/*                  ATAN                        */
/************************************************/
template <class T>
canonical_polynomial<T> atan(const canonical_polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    canonical_polynomial<T> res(nvar,degree);

    //Canonical expansion of atan in [a,b]
    std::vector<T> range = other.get_range();

    //approximate atan [-range,range]
    std::vector<T> cheb_atan = canonical_polynomial<T>::approximation_1d(atan,range[0],range[1],other.get_degree());
    //univariate composition
    std::vector<canonical_polynomial<T> > base = other.evaluate_base(range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_atan[i];
    }

    return res;
}
template class canonical_polynomial<double>
atan(const canonical_polynomial<double> &);
template class canonical_polynomial<float>
atan(const canonical_polynomial<float> &);
template class canonical_polynomial<long double>
atan(const canonical_polynomial<long double> &);

/************************************************/
/*                  ACOT                        */
/************************************************/
template <class T>
canonical_polynomial<T> acot(const canonical_polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;
    return canonical_polynomial<T>(0,0);
}
template class canonical_polynomial<double>
acot(const canonical_polynomial<double> &);
template class canonical_polynomial<float>
acot(const canonical_polynomial<float> &);
template class canonical_polynomial<long double>
acot(const canonical_polynomial<long double> &);



// OTHERS
//EXPONENTIAL FUNCTION
template <class T>
/************************************************/
/*                  EXP                         */
/************************************************/
canonical_polynomial<T> exp(const canonical_polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    canonical_polynomial<T> res(nvar,degree);

    //Canonical expansion of exp in [a,b]
    std::vector<T> range = other.get_range();

    //approximate exp [-range,range]
    std::vector<T> cheb_exp = canonical_polynomial<T>::approximation_1d(exp,range[0],range[1],other.get_degree());
    //univariate composition
    std::vector<canonical_polynomial<T> > base = other.evaluate_base(range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_exp[i];
    }

    return res;
}
template class canonical_polynomial<double>
exp(const canonical_polynomial<double> &);
template class canonical_polynomial<float>
exp(const canonical_polynomial<float> &);
template class canonical_polynomial<long double>
exp(const canonical_polynomial<long double> &);


/************************************************/
/*                  SQRT                        */
/************************************************/
template <class T>
canonical_polynomial<T> sqrt(const canonical_polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    canonical_polynomial<T> res(nvar,degree);

    //Canonical expansion of exp in [a,b]
    std::vector<T> range = other.get_range();

    //approximate exp [-range,range]
    std::vector<T> cheb_sqrt = canonical_polynomial<T>::approximation_1d(sqrt,range[0],range[1],other.get_degree());
    //univariate composition
    std::vector<canonical_polynomial<T> > base = other.evaluate_base(range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_sqrt[i];
    }

    return res;
}
template class canonical_polynomial<double>
sqrt(const canonical_polynomial<double> &);
template class canonical_polynomial<float>
sqrt(const canonical_polynomial<float> &);
template class canonical_polynomial<long double>
sqrt(const canonical_polynomial<long double> &);

/************************************************/
/*                  LOG                         */
/************************************************/
template <class T>
canonical_polynomial<T> log(const canonical_polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    canonical_polynomial<T> res(nvar,degree);

    //Canonical expansion of exp in [a,b]
    std::vector<T> range = other.get_range();

    //approximate exp [-range,range]
    std::vector<T> cheb_log = canonical_polynomial<T>::approximation_1d(log,range[0],range[1],other.get_degree());
    //univariate composition
    std::vector<canonical_polynomial<T> > base = other.evaluate_base(range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_log[i];
    }

    return res;
}
template class canonical_polynomial<double>
log(const canonical_polynomial<double> &);
template class canonical_polynomial<float>
log(const canonical_polynomial<float> &);
template class canonical_polynomial<long double>
log(const canonical_polynomial<long double> &);

/************************************************/
/*                  LOG10                       */
/************************************************/
template <class T>
canonical_polynomial<T> log10(const canonical_polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;
        return canonical_polynomial<T>(0,0);
}
template class canonical_polynomial<double>
log10(const canonical_polynomial<double> &);
template class canonical_polynomial<float>
log10(const canonical_polynomial<float> &);
template class canonical_polynomial<long double>
log10(const canonical_polynomial<long double> &);

/************************************************/
/*                  POW                         */
/************************************************/
template <class T>
canonical_polynomial<T> pow(const canonical_polynomial<T> &other, const int &exponent){
    if(exponent<=1){
        std::cout<<"pow function with integer exponent, integer must be > 1"<<std::endl;
        exit(EXIT_FAILURE);
    }
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    canonical_polynomial<T> res(nvar,degree);

    res = other;
    for(int i=1; i<exponent; i++)
        res *= other;

    return res;
}
template class canonical_polynomial<double>
pow(const canonical_polynomial<double> &, const int &);
template class canonical_polynomial<float>
pow(const canonical_polynomial<float> &, const int &);
template class canonical_polynomial<long double>
pow(const canonical_polynomial<long double> &, const int &);



template <class T>
canonical_polynomial<T> pow(const canonical_polynomial<T> &other, const double &exponent){
std::cout<<"NOT IMPLEMENTED"<<std::endl;
    return canonical_polynomial<T>(0,0);
}
template class canonical_polynomial<double>
pow(const canonical_polynomial<double> &, const double &);
template class canonical_polynomial<float>
pow(const canonical_polynomial<float> &, const double &);
template class canonical_polynomial<long double>
pow(const canonical_polynomial<long double> &, const double &);
