#include "Intrusive/elementary_functions.h"

using namespace smart;
using namespace intrusive;

/************************************************/
/*           DIRECT MULTIPLICATION              */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> direct_multiplication(const Chebyshev_Polynomial<T> &x0, const Chebyshev_Polynomial<T> &x1) {
    if(x0.get_nvar()!=x1.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(x0.get_degree()!=x1.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    Chebyshev_Polynomial<T> res(x0.get_nvar(),x0.get_degree());
    std::vector<T> res_coeffs(combination(x0.get_nvar(),x0.get_degree()));
    double nvariations = pow(2,x0.get_nvar());
    std::vector<T> x0_coeffs = x0.get_coeffs();
    std::vector<T> x1_coeffs = x1.get_coeffs();
    std::vector<std::vector<int> > x0_J=x0.get_J();
    std::vector<std::vector<int> > x0_N=x0.get_N();
    std::vector<std::vector<int> > x0_t=x0.get_t();
    for(int i=0; i<=x0.get_degree(); i++){//loop over subset degree i of poly1
        for(int j=0; j<=x1.get_degree(); j++){//loop over subset degree j of poly2
            //if((i+j)<=m_degree){
                for(int idx1=0; idx1<x0_J[x0.get_nvar()][i]; idx1++){//index over elements with degree i in poly1
                    for(int idx2=0; idx2<x0_J[x0.get_nvar()][j]; idx2++){//index over elements with degree j in poly2
                        int sub_idx1=0, sub_idx2=0, sub_idx3=0;
                        if(i>0) sub_idx2=x0_N[x0.get_nvar()][i-1];
                        if(j>0) sub_idx3=x0_N[x0.get_nvar()][j-1];
                        if(fabs(x0_coeffs[sub_idx2+idx1])>ZERO && fabs(x1_coeffs[sub_idx3+idx2])>ZERO){
                            std::vector<int> v1 = x0.get_row(idx1,i);
                            std::vector<int> v2 = x0.get_row(idx2,j);
                            std::vector<int> v3(x0.get_nvar());
                            for(int iter=0; iter<nvariations; iter++){
                                for(int k=0; k<x0.get_nvar(); k++){
                                    v3[k] = std::fabs(v1[k]+x0_t[iter][k]*v2[k]);
                                }
                                int deg3 = std::accumulate(v3.begin(),v3.end(),0);
                                if(deg3<=x0.get_degree()){
                                    int pos = res.get_idx(v3);
                                    sub_idx1 = 0;
                                    if(deg3>0) sub_idx1=x0_N[x0.get_nvar()][deg3-1];
                                    res_coeffs[sub_idx1 + pos] +=
                                        (1.0/nvariations)*(x0_coeffs[sub_idx2+idx1]*x1_coeffs[sub_idx3+idx2]);
                                }
                            }
                        }
                    }
                }
            //}
        }
    }

    res.set_coeffs(res_coeffs);
    return res;
}
template class Chebyshev_Polynomial<double>
direct_multiplication(const Chebyshev_Polynomial<double> &, const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
direct_multiplication(const Chebyshev_Polynomial<float> &, const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
direct_multiplication(const Chebyshev_Polynomial<long double> &, const Chebyshev_Polynomial<long double> &);


/************************************************/
/*                  SIN                         */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> sin(const Chebyshev_Polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Chebyshev_Polynomial<T> res(nvar,degree);

    //chebyshev expansion of sin in [a,b]
    std::vector<T> range = other.get_range();

    //approximate sin [-range,range]
    std::vector<T> cheb_sin = Chebyshev_Polynomial<T>::cheb_approximation(sin,range[0],range[1]);
    //univariate composition
    std::vector<Chebyshev_Polynomial<T> > base = Chebyshev_Polynomial<T>::evaluate_base(other,range[0],range[1]);
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
    std::vector<T> range = other.get_range();

    //approximate cos [-range,range]
    std::vector<T> cheb_cos = Chebyshev_Polynomial<T>::cheb_approximation(cos,range[0],range[1]);
    //univariate composition
    std::vector<Chebyshev_Polynomial<T> > base = Chebyshev_Polynomial<T>::evaluate_base(other,range[0],range[1]);
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
    std::vector<T> range = other.get_range();

    //approximate tan [-range,range]
    std::vector<T> cheb_tan = Chebyshev_Polynomial<T>::cheb_approximation(tan,range[0],range[1]);
    //univariate composition
    std::vector<Chebyshev_Polynomial<T> > base = Chebyshev_Polynomial<T>::evaluate_base(other,range[0],range[1]);
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
    std::vector<T> range = other.get_range();

    //approximate atan [-range,range]
    std::vector<T> cheb_atan = Chebyshev_Polynomial<T>::cheb_approximation(atan,range[0],range[1]);
    //univariate composition
    std::vector<Chebyshev_Polynomial<T> > base = Chebyshev_Polynomial<T>::evaluate_base(other,range[0],range[1]);
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
    std::vector<T> range = other.get_range();

    //approximate exp [-range,range]
    std::vector<T> cheb_exp = Chebyshev_Polynomial<T>::cheb_approximation(exp,range[0],range[1]);
    //univariate composition
    std::vector<Chebyshev_Polynomial<T> > base = Chebyshev_Polynomial<T>::evaluate_base(other,range[0],range[1]);
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
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Chebyshev_Polynomial<T> res(nvar,degree);

    //chebyshev expansion of exp in [a,b]
    std::vector<T> range = other.get_range();

    //approximate exp [-range,range]
    std::vector<T> cheb_sqrt = Chebyshev_Polynomial<T>::cheb_approximation(sqrt,range[0],range[1]);
    //univariate composition
    std::vector<Chebyshev_Polynomial<T> > base = Chebyshev_Polynomial<T>::evaluate_base(other,range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_sqrt[i];
    }

    return res;
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
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Chebyshev_Polynomial<T> res(nvar,degree);

    //chebyshev expansion of exp in [a,b]
    std::vector<T> range = other.get_range();

    //approximate exp [-range,range]
    std::vector<T> cheb_log = Chebyshev_Polynomial<T>::cheb_approximation(log,range[0],range[1]);
    //univariate composition
    std::vector<Chebyshev_Polynomial<T> > base = Chebyshev_Polynomial<T>::evaluate_base(other,range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_log[i];
    }

    return res;
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