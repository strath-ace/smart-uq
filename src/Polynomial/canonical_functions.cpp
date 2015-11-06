#include "Polynomial/canonical_functions.h"

using namespace smart;
using namespace polynomial;

/************************************************/
/*                  SIN                         */
/************************************************/
template <class T>
Canonical_Polynomial<T> sin(const Canonical_Polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Canonical_Polynomial<T> res(nvar,degree);

    //Canonical expansion of sin in [a,b]
    std::vector<T> range = other.get_range();

    //approximate sin [-range,range]
    std::vector<T> cheb_sin = Canonical_Polynomial<T>::approximation_1d(sin,range[0],range[1],other.get_degree());
    //univariate composition
    std::vector<Canonical_Polynomial<T> > base = other.evaluate_base(range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_sin[i];
    }

    return res;
}
template class Canonical_Polynomial<double>
sin(const Canonical_Polynomial<double> &);
template class Canonical_Polynomial<float>
sin(const Canonical_Polynomial<float> &);
template class Canonical_Polynomial<long double>
sin(const Canonical_Polynomial<long double> &);


/************************************************/
/*                  COS                         */
/************************************************/
template <class T>
Canonical_Polynomial<T> cos(const Canonical_Polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Canonical_Polynomial<T> res(nvar,degree);

    //Canonical expansion of cos in [a,b]
    std::vector<T> range = other.get_range();

    //approximate cos [-range,range]
    std::vector<T> cheb_cos = Canonical_Polynomial<T>::approximation_1d(cos,range[0],range[1],other.get_degree());
    //univariate composition
    std::vector<Canonical_Polynomial<T> > base = other.evaluate_base(range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_cos[i];
    }

    return res;
}
template class Canonical_Polynomial<double>
cos(const Canonical_Polynomial<double> &);
template class Canonical_Polynomial<float>
cos(const Canonical_Polynomial<float> &);
template class Canonical_Polynomial<long double>
cos(const Canonical_Polynomial<long double> &);

/************************************************/
/*                  TAN                         */
/************************************************/
template <class T>
Canonical_Polynomial<T> tan(const Canonical_Polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Canonical_Polynomial<T> res(nvar,degree);

    //Canonical expansion of tan in [a,b]
    std::vector<T> range = other.get_range();

    //approximate tan [-range,range]
    std::vector<T> cheb_tan = Canonical_Polynomial<T>::approximation_1d(tan,range[0],range[1],other.get_degree());
    //univariate composition
    std::vector<Canonical_Polynomial<T> > base = other.evaluate_base(range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_tan[i];
    }

    return res;
}
template class Canonical_Polynomial<double>
tan(const Canonical_Polynomial<double> &);
template class Canonical_Polynomial<float>
tan(const Canonical_Polynomial<float> &);
template class Canonical_Polynomial<long double>
tan(const Canonical_Polynomial<long double> &);

/************************************************/
/*                  COT                         */
/************************************************/
template <class T>
Canonical_Polynomial<T> cot(const Canonical_Polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;

}
template class Canonical_Polynomial<double>
cot(const Canonical_Polynomial<double> &);
template class Canonical_Polynomial<float>
cot(const Canonical_Polynomial<float> &);
template class Canonical_Polynomial<long double>
cot(const Canonical_Polynomial<long double> &);

/************************************************/
/*                  ASIN                        */
/************************************************/
template <class T>
Canonical_Polynomial<T> asin(const Canonical_Polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;

}
template class Canonical_Polynomial<double>
asin(const Canonical_Polynomial<double> &);
template class Canonical_Polynomial<float>
asin(const Canonical_Polynomial<float> &);
template class Canonical_Polynomial<long double>
asin(const Canonical_Polynomial<long double> &);

/************************************************/
/*                  ACOS                        */
/************************************************/
template <class T>
Canonical_Polynomial<T> acos(const Canonical_Polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;

}
template class Canonical_Polynomial<double>
acos(const Canonical_Polynomial<double> &);
template class Canonical_Polynomial<float>
acos(const Canonical_Polynomial<float> &);
template class Canonical_Polynomial<long double>
acos(const Canonical_Polynomial<long double> &);

/************************************************/
/*                  ATAN                        */
/************************************************/
template <class T>
Canonical_Polynomial<T> atan(const Canonical_Polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Canonical_Polynomial<T> res(nvar,degree);

    //Canonical expansion of atan in [a,b]
    std::vector<T> range = other.get_range();

    //approximate atan [-range,range]
    std::vector<T> cheb_atan = Canonical_Polynomial<T>::approximation_1d(atan,range[0],range[1],other.get_degree());
    //univariate composition
    std::vector<Canonical_Polynomial<T> > base = other.evaluate_base(range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_atan[i];
    }

    return res;
}
template class Canonical_Polynomial<double>
atan(const Canonical_Polynomial<double> &);
template class Canonical_Polynomial<float>
atan(const Canonical_Polynomial<float> &);
template class Canonical_Polynomial<long double>
atan(const Canonical_Polynomial<long double> &);

/************************************************/
/*                  ACOT                        */
/************************************************/
template <class T>
Canonical_Polynomial<T> acot(const Canonical_Polynomial<T> &other){
    std::cout<<"NOT IMPLEMENTED"<<std::endl;

}
template class Canonical_Polynomial<double>
acot(const Canonical_Polynomial<double> &);
template class Canonical_Polynomial<float>
acot(const Canonical_Polynomial<float> &);
template class Canonical_Polynomial<long double>
acot(const Canonical_Polynomial<long double> &);



// OTHERS
//EXPONENTIAL FUNCTION
template <class T>
/************************************************/
/*                  EXP                         */
/************************************************/
Canonical_Polynomial<T> exp(const Canonical_Polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Canonical_Polynomial<T> res(nvar,degree);

    //Canonical expansion of exp in [a,b]
    std::vector<T> range = other.get_range();

    //approximate exp [-range,range]
    std::vector<T> cheb_exp = Canonical_Polynomial<T>::approximation_1d(exp,range[0],range[1],other.get_degree());
    //univariate composition
    std::vector<Canonical_Polynomial<T> > base = other.evaluate_base(range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_exp[i];
    }

    return res;
}
template class Canonical_Polynomial<double>
exp(const Canonical_Polynomial<double> &);
template class Canonical_Polynomial<float>
exp(const Canonical_Polynomial<float> &);
template class Canonical_Polynomial<long double>
exp(const Canonical_Polynomial<long double> &);


/************************************************/
/*                  SQRT                        */
/************************************************/
template <class T>
Canonical_Polynomial<T> sqrt(const Canonical_Polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Canonical_Polynomial<T> res(nvar,degree);

    //Canonical expansion of exp in [a,b]
    std::vector<T> range = other.get_range();

    //approximate exp [-range,range]
    std::vector<T> cheb_sqrt = Canonical_Polynomial<T>::approximation_1d(sqrt,range[0],range[1],other.get_degree());
    //univariate composition
    std::vector<Canonical_Polynomial<T> > base = other.evaluate_base(range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_sqrt[i];
    }

    return res;
}
template class Canonical_Polynomial<double>
sqrt(const Canonical_Polynomial<double> &);
template class Canonical_Polynomial<float>
sqrt(const Canonical_Polynomial<float> &);
template class Canonical_Polynomial<long double>
sqrt(const Canonical_Polynomial<long double> &);

/************************************************/
/*                  LOG                         */
/************************************************/
template <class T>
Canonical_Polynomial<T> log(const Canonical_Polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Canonical_Polynomial<T> res(nvar,degree);

    //Canonical expansion of exp in [a,b]
    std::vector<T> range = other.get_range();

    //approximate exp [-range,range]
    std::vector<T> cheb_log = Canonical_Polynomial<T>::approximation_1d(log,range[0],range[1],other.get_degree());
    //univariate composition
    std::vector<Canonical_Polynomial<T> > base = other.evaluate_base(range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_log[i];
    }

    return res;
}
template class Canonical_Polynomial<double>
log(const Canonical_Polynomial<double> &);
template class Canonical_Polynomial<float>
log(const Canonical_Polynomial<float> &);
template class Canonical_Polynomial<long double>
log(const Canonical_Polynomial<long double> &);

/************************************************/
/*                  LOG10                       */
/************************************************/
template <class T>
Canonical_Polynomial<T> log10(const Canonical_Polynomial<T> &other){

}
template class Canonical_Polynomial<double>
log10(const Canonical_Polynomial<double> &);
template class Canonical_Polynomial<float>
log10(const Canonical_Polynomial<float> &);
template class Canonical_Polynomial<long double>
log10(const Canonical_Polynomial<long double> &);

/************************************************/
/*                  POW                         */
/************************************************/
template <class T>
Canonical_Polynomial<T> pow(const Canonical_Polynomial<T> &other, const int &exponent){
    if(exponent<=1){
        std::cout<<"pow function with integer exponent, integer must be > 1"<<std::endl;
        exit(EXIT_FAILURE);
    }
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Canonical_Polynomial<T> res(nvar,degree);

    res = other;
    for(int i=1; i<exponent; i++)
        res *= other;

    return res;
}
template class Canonical_Polynomial<double>
pow(const Canonical_Polynomial<double> &, const int &);
template class Canonical_Polynomial<float>
pow(const Canonical_Polynomial<float> &, const int &);
template class Canonical_Polynomial<long double>
pow(const Canonical_Polynomial<long double> &, const int &);



template <class T>
Canonical_Polynomial<T> pow(const Canonical_Polynomial<T> &other, const double &exponent){
std::cout<<"NOT IMPLEMENTED"<<std::endl;
}
template class Canonical_Polynomial<double>
pow(const Canonical_Polynomial<double> &, const double &);
template class Canonical_Polynomial<float>
pow(const Canonical_Polynomial<float> &, const double &);
template class Canonical_Polynomial<long double>
pow(const Canonical_Polynomial<long double> &, const double &);