#include "elementary_functions.h"

template <class T>
Chebyshev_Polynomial<T> composition(const std::vector<T> &coeffs, const Chebyshev_Polynomial<T> &other, const T &range){

    int nvar =  other.get_nvar();
    int degree = other.get_degree();

    std::vector<Chebyshev_Polynomial<T> > v;
    Chebyshev_Polynomial<T> res(nvar, degree);

    for(int i=0; i<=degree; i++){
        v.push_back(Chebyshev_Polynomial<T>(nvar,degree));
    }

    v[0] = 1.0;
    v[1] = other/range;

    for (int i=2; i<=degree; i++){
        v[i] = 2.0 * other/range * v[i-1] - v[i-2];
    }

    for (int i=0; i<=degree; i++){
        res += v[i]*coeffs[i];
    }

    return res;
}

template <class T>
std::vector<T> cheb_approximation(T (*f)(T x), const T a, const T b){
    int n = Chebyshev_Polynomial<T>::MAX_DEGREE;
    std::vector<T> res(n+1), d(n+1);
    T fac;
    T pi = 3.141592653589793;
    T t;
    T total;
    T y;

    for (int k = 0; k <= n; k++)
    {
      t = cos(pi*(k+0.5)/(n+1)); //zeros Ti
      y = ((1.0+t)*b + (1.0-t)*a)/2.0; //mapped zeros
      d[k] = f(y); //evaluate function
    }

    fac = 2.0/(n+1);

    for (int j = 0; j <= n; j++)
    {
      total = 0.0;
      for (int k = 0; k <= n; k++)
      {
        total = total+d[k]*cos( (pi*j)*( (k+ 0.5)/(n+1) ) );
      }
      res[j] = fac*total;
    }

    res[0] = res[0]/2.0;

    return res;
}
template class std::vector<double>
cheb_approximation(double (*f)(double x), const double a, const double b);
template class std::vector<float>
cheb_approximation(float (*f)(float x), const float a, const float b);
template class std::vector<long double>
cheb_approximation(long double (*f)(long double x), const long double a, const long double b);

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
    std::vector<T> cheb_sin = cheb_approximation(sin,-range,range);

    res = composition(cheb_sin, other,range);

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
    T a = -1.0;
    T b = 1.0;
    std::vector<T> cheb_cos = cheb_approximation(cos,a,b);

    res = composition(cheb_cos, other);

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
    T a = -1.0;
    T b = 1.0;
    std::vector<T> cheb_tan = cheb_approximation(tan,a,b);

    res = composition(cheb_tan, other);

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
    T a = -1.0;
    T b = 1.0;
    std::vector<T> cheb_atan = cheb_approximation(atan,a,b);

    res = composition(cheb_atan, other);

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
    T a = -1.0;
    T b = 1.0;
    std::vector<T> cheb_exp = cheb_approximation(exp,a,b);

    res = composition(cheb_exp, other);

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


/************************************************/
/*                  INV                         */
/************************************************/
template <class T>
Chebyshev_Polynomial<T> inv(const Chebyshev_Polynomial<T> &other){
std::cout<<"NOT IMPLEMENTED"<<std::endl;

}
template class Chebyshev_Polynomial<double>
inv(const Chebyshev_Polynomial<double> &);
template class Chebyshev_Polynomial<float>
inv(const Chebyshev_Polynomial<float> &);
template class Chebyshev_Polynomial<long double>
inv(const Chebyshev_Polynomial<long double> &);
