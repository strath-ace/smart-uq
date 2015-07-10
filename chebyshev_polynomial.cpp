/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
---------------- e-mail: annalisa.riccardi@strath.ac.uk ----------------------
------------------------- Author: Annalisa Riccardi --------------------------
*/


#include "chebyshev_polynomial.h"

template <class T>
Chebyshev_Polynomial<T>::Chebyshev_Polynomial(const int &nvar, const int &order): m_coeffs(0), m_degree(0), m_nvar(0){
    //allocate memory for coefficients vector

    if(order > MAX_DEGREE){
        std::cout<<"Maximum allowed polynomial degree is 20";
        exit(EXIT_FAILURE);
    }

    int n = combination(nvar,order);
    m_coeffs.resize(n);

    //save some info
    m_degree = order;
    m_nvar = nvar;

    m_J.resize(nvar+1);
    m_N.resize(nvar+1);
    m_t.resize(pow(2,nvar));
    for(int i=0; i<=nvar; i++){
        m_J[i].resize(order+1);
        m_N[i].resize(order+1);
        if(i<nvar)
            m_t[i].resize(2);
    }

    initialize_J();
    initialize_N();
    initialize_t();

}

template <class T>
Chebyshev_Polynomial<T>::Chebyshev_Polynomial(const int &nvar, const int &order, const int &i): m_coeffs(0), m_degree(0), m_nvar(0){

    if(order > MAX_DEGREE){
        std::cout<<"Maximum allowed polynomial degree is 20";
        exit(EXIT_FAILURE);
    }

    //save some info
    m_degree = order;
    m_nvar = nvar;

    if(i>=m_nvar){
        std::cout<<"base elements index are from [0,nvar-1]";
        exit(EXIT_FAILURE);
    }
    else{
        //allocate memory for coefficients vector

        int n = combination(nvar,order);
        m_coeffs.resize(n);
        m_coeffs[i+1] = 1.0;

        m_J.resize(nvar+1);
        m_N.resize(nvar+1);
        m_t.resize(pow(2,nvar));
        for(int i=0; i<=nvar; i++){
            m_J[i].resize(order+1);
            m_N[i].resize(order+1);
            if(i<nvar)
                m_t[i].resize(2);
        }

        initialize_J();
        initialize_N();
        initialize_t();

    }

}

template <class T>
Chebyshev_Polynomial<T>::Chebyshev_Polynomial(const int &nvar, const int &order, const T &value): m_coeffs(0), m_degree(0), m_nvar(0){

    if(order > MAX_DEGREE){
        std::cout<<"Maximum allowed polynomial degree is 20";
        exit(EXIT_FAILURE);
    }

    //save some info
    m_degree = order;
    m_nvar = nvar;

    //allocate memory for coefficients vector

    int n = combination(nvar,order);
    m_coeffs.resize(n);
    m_coeffs[0] = value;

    m_J.resize(nvar+1);
    m_N.resize(nvar+1);
    m_t.resize(pow(2,nvar));
    for(int i=0; i<=nvar; i++){
        m_J[i].resize(order+1);
        m_N[i].resize(order+1);
        if(i<nvar)
            m_t[i].resize(2);
    }

    initialize_J();
    initialize_N();
    initialize_t();


}


template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator+(const Chebyshev_Polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(m_degree!=other.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    int n = get_coeffs().size();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    Chebyshev_Polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs[i] + other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator-(const Chebyshev_Polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(m_degree!=other.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    int n = get_coeffs().size();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    Chebyshev_Polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs[i] - other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator*(const Chebyshev_Polynomial<T> &other) const{
    if(m_nvar!=other.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(m_degree!=other.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    Chebyshev_Polynomial<T> res(m_nvar,m_degree);
    std::vector<T> res_coeffs(combination(m_nvar,m_degree));
    double nvariations = pow(2,m_nvar);
    std::vector<T> other_coeffs = other.get_coeffs();
    for(int i=0; i<=m_degree; i++){//loop over subset degree i of poly1
        for(int j=0; j<=other.get_degree(); j++){//loop over subset degree j of poly2
            //if((i+j)<=m_degree){
                for(int idx1=0; idx1<m_J[m_nvar][i]; idx1++){//index over elements with degree i in poly1
                    for(int idx2=0; idx2<m_J[m_nvar][j]; idx2++){//index over elements with degree j in poly2
                        int sub_idx1=0, sub_idx2=0, sub_idx3=0;
                        if(i>0) sub_idx2=m_N[m_nvar][i-1];
                        if(j>0) sub_idx3=m_N[m_nvar][j-1];
                        if(fabs(m_coeffs[sub_idx2+idx1])>ZERO && fabs(other_coeffs[sub_idx3+idx2])>ZERO){
                            std::vector<int> v1 = get_row(idx1,i);
                            std::vector<int> v2 = get_row(idx2,j);
                            std::vector<int> v3(m_nvar);
                            for(int iter=0; iter<nvariations; iter++){
                                for(int k=0; k<m_nvar; k++){
                                    v3[k] = std::fabs(v1[k]+m_t[iter][k]*v2[k]);
                                }
                                int deg3 = std::accumulate(v3.begin(),v3.end(),0);
                                if(deg3<=m_degree){
                                    int pos = res.get_idx(v3);
                                    sub_idx1 = 0;
                                    if(deg3>0) sub_idx1=m_N[m_nvar][deg3-1];
                                    res_coeffs[sub_idx1 + pos] +=
                                        (1.0/nvariations)*(m_coeffs[sub_idx2+idx1]*other_coeffs[sub_idx3+idx2]);
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

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::inv(const Chebyshev_Polynomial<T> &other) const{
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    Chebyshev_Polynomial<T> res(nvar,degree);

    //chebyshev expansion of sin in [a,b]
    //    std::vector<T> coeffs = other.get_coeffs();
    std::vector<T> range = other.get_range();
    T a, b;

    //    //computing value of other in zero
    //    T value = 0.0;
    //    int count = 0;
    //    for(int deg=0; deg<=degree; deg++){
    //        for(int i=0; i<m_J[nvar][deg]; i++){
    //            std::vector<int> row = get_row(i,deg); //get for example vector (1 0 0) = x, (0 1 0) = y...
    //            int prod = 1;
    //            for(int j=0;j<m_nvar; j++){
    //                if(row[j]%2 == 0){
    //                    if(row[j]%4 != 0)
    //                        prod*= -1.0;
    //                }
    //                else prod *= 0.0;
    //            }
    //            value += coeffs[count]*prod;
    //            count++;
    //        }
    //    }

    //    T one = 1.0;
    //    if(value>0){
    //        a = min(one,range[1]);
    //        b = max(one,range[1]);
    //    }
    //    else if(value<0){
    //        a = max(-one,range[0]);
    //        b = min(-one,range[0]);
    //    }
    //    else{
    //        std::cout<<"Inverting a function with a zero in the interval [-1,1]."<<std::endl;
    //        exit(EXIT_FAILURE);
    //    }

    a = range[0];
    b = range[1];
    std::vector<T> cheb_inv = cheb_approximation(inverse,a,b);
    //univariate composition
    std::vector<Chebyshev_Polynomial<T> > base = evaluate_base(other, a,b);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_inv[i];
    }

    return res;
}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator/(const Chebyshev_Polynomial<T> &other) const{
    if(m_nvar!=other.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(m_degree!=other.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    Chebyshev_Polynomial<T> res(m_nvar,m_degree);
    res = inv(other);

    return res*(*this);
}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator+(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    Chebyshev_Polynomial<T> res(m_nvar,m_degree);

    coeffs[0] += other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator-(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    Chebyshev_Polynomial<T> res(m_nvar,m_degree);

    coeffs[0] -= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator*(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    Chebyshev_Polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<coeffs.size(); i++)
        coeffs[i] *= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator/(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    Chebyshev_Polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<coeffs.size(); i++)
        coeffs[i] /= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Chebyshev_Polynomial<T>& Chebyshev_Polynomial<T>::operator=(const Chebyshev_Polynomial<T> &other){
    if(m_nvar!=other.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(m_degree!=other.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    m_coeffs = other.get_coeffs();
    return *this;
}

template <class T>
Chebyshev_Polynomial<T>& Chebyshev_Polynomial<T>::operator=(const T &other){

    std::vector<T> coeffs(m_coeffs.size());
    coeffs[0] = other;

    m_coeffs = coeffs;
    return *this;
}

template <class T>
Chebyshev_Polynomial<T>& Chebyshev_Polynomial<T>::operator+=(const Chebyshev_Polynomial<T> &other){
    *this = Chebyshev_Polynomial<T>::operator+(other);
    return *this;
}

template <class T>
Chebyshev_Polynomial<T>& Chebyshev_Polynomial<T>::operator-=(const Chebyshev_Polynomial<T> &other){
    *this = Chebyshev_Polynomial<T>::operator-(other);
    return *this;
}

template <class T>	
Chebyshev_Polynomial<T>& Chebyshev_Polynomial<T>::operator*=(const Chebyshev_Polynomial<T> &other){
    *this = Chebyshev_Polynomial<T>::operator*(other);
    return *this;
}

template <class T>
Chebyshev_Polynomial<T>& Chebyshev_Polynomial<T>::operator/=(const Chebyshev_Polynomial<T> &other){
    *this = Chebyshev_Polynomial<T>::operator/(other);
    return *this;
}

template <class T>
Chebyshev_Polynomial<T>& Chebyshev_Polynomial<T>::operator+=(const T& other){
    *this = Chebyshev_Polynomial<T>::operator+(other);
    return *this;
}

template <class T>
Chebyshev_Polynomial<T>& Chebyshev_Polynomial<T>::operator-=(const T& other){
    *this = Chebyshev_Polynomial<T>::operator-(other);
    return *this;
}

template <class T>
Chebyshev_Polynomial<T>& Chebyshev_Polynomial<T>::operator*=(const T& other){
    *this = Chebyshev_Polynomial<T>::operator*(other);
    return *this;
}

template <class T>
Chebyshev_Polynomial<T>& Chebyshev_Polynomial<T>::operator/=(const T& other){
    *this = Chebyshev_Polynomial<T>::operator/(other);
    return *this;
}

template <class T>
bool Chebyshev_Polynomial<T>::operator==(const Chebyshev_Polynomial<T> &other) const{
    if(m_nvar!=other.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(m_degree!=other.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    if(m_coeffs==other.get_coeffs())
        return true;
    return false;
}

template <class T>	
bool Chebyshev_Polynomial<T>::operator!=(const Chebyshev_Polynomial<T> &other) const{
    if(m_nvar!=other.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(m_degree!=other.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    if(Chebyshev_Polynomial<T>::operator==(other)) return false;
    else return true;

}


//evaluate chebyshev base t0(x), t1(x), t2(x) in a polynomial. It first map x from [a,b] to [-1,1]
template <class T>
std::vector<Chebyshev_Polynomial<T> > Chebyshev_Polynomial<T>::evaluate_base(const Chebyshev_Polynomial<T> &other, const T &a, const T &b){

    int nvar =  other.get_nvar();
    int degree = other.get_degree();

    std::vector<Chebyshev_Polynomial<T> > v;

    for(int i=0; i<=degree; i++){
        v.push_back(Chebyshev_Polynomial<T>(nvar,degree));
    }

    //mapping the argument from the domain of composition to [-1,1]
    Chebyshev_Polynomial<T> mapped(nvar,degree);
    if(b==a)
        mapped.set_coeffs(0,a);
    else
        mapped = (2.0*other-(a+b))/(b-a);

    v[0] = 1.0;
    v[1] = mapped;

    for (int i=2; i<=degree; i++){
        v[i] = 2.0 * mapped * v[i-1] - v[i-2];
    }

    return v;
}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::composition(const std::vector<Chebyshev_Polynomial<T> > &other) const{
    if(m_nvar!=other.size()){
        std::cout<<"Composition is with a vector of polynomial of the same size of nvar"<<std::endl;
        exit(EXIT_FAILURE);
    }
    for(int i=0; i<m_nvar; i++){
        if(m_nvar!=other[i].get_nvar()){
            std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
            exit(EXIT_FAILURE);
        }
        if(m_degree!=other[i].get_degree()){
            std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
            exit(EXIT_FAILURE);
        }
    }

    //allocate memory
    std::vector<std::vector<Chebyshev_Polynomial<T> > > base;
    for(int j=0; j<m_nvar; j++){
        std::vector<Chebyshev_Polynomial<T> > v;
        for(int i=0; i<=m_degree; i++){
            v.push_back(Chebyshev_Polynomial<T>(m_nvar,m_degree));
        }
        base.push_back(v);
    }

    //evaluate all basis
    for(int j=0; j<m_nvar; j++){
        //T range = other[j].get_range();
        std::vector<Chebyshev_Polynomial<T> > v = evaluate_base(other[j],-1.0,1.0);
        for(int i=0; i<=m_degree; i++)
            base[j][i] = v[i];
    }

    //composing
    Chebyshev_Polynomial<T> res(m_nvar,m_degree);
    int count = 0;
    for(int deg=0; deg<=m_degree; deg++){
        for(int i=0; i<m_J[m_nvar][deg]; i++){
            std::vector<int> row = get_row(i,deg); //get for example vector (1 0 0) = x, (0 1 0) = y...
            Chebyshev_Polynomial<T> prod(m_nvar,m_degree);
            prod.set_coeffs(0,1.0);
            for(int j=0;j<m_nvar; j++){
                prod*=base[j][row[j]];
            }
            res += m_coeffs[count]*prod;
            count++;
        }
    }

    return res;
}

template <class T>
std::vector<T> Chebyshev_Polynomial<T>::cheb_approximation(T (*f)(T x), const T a, const T b){
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


//not part of the algebra, private routines
template <class T>
void Chebyshev_Polynomial<T>::initialize_J()
{
    int i,j;

    //fill J
    for(j = 0; j <= m_degree; ++j)
        m_J[1][j] = 1;
    for(i = 2; i <= m_nvar; ++i) {
        m_J[i][0] = 1;
        for(j = 1; j <= m_degree; ++j)
            m_J[i][j] = m_J[i][j-1]+m_J[i-1][j];
    }

}

template <class T>
void Chebyshev_Polynomial<T>::initialize_N()
{
    int i,j;
    //fill N
    for(i = 1; i <= m_nvar; ++i)
        m_N[i][0] = m_J[i][0];
    for(i = 1; i <= m_nvar; ++i) {
        for(j = 1; j <= m_degree; ++j)
            m_N[i][j] = m_N[i][j-1]+m_J[i][j];
    }

}

template <class T>
void Chebyshev_Polynomial<T>::initialize_t(){
    std::vector<int> values(2);
    values[0] = -1;
    values[1] = 1;

    variations(values,m_nvar, m_t);

}

template <class T>
int Chebyshev_Polynomial<T>::get_idx(const std::vector<int> &k) const
{
    int i, j, l;

    int idx = 0;
    for(i = m_nvar-1; (i>=0) && (k[i]==0); --i);

    if(i < 0)
        return idx;
    j = m_nvar-i;
    l = -1;
    while(i >= 0) {
        l += k[i--];
        idx += m_N[j++][l];
    }
    idx -= m_N[m_nvar][l];

    return idx;
}

template <class T>
std::vector<int> Chebyshev_Polynomial<T>::get_row(const int &idx, const int &deg) const
{
    int i, j, l;
    int idx_tmp = idx;
    std::vector<int> k(m_nvar);

    if(deg > 0)
        idx_tmp += m_N[m_nvar][deg-1];
    l = m_nvar;
    for(i = 0; i < m_nvar; ++i) {
        if(idx_tmp == 0) {
            for(j = i; j < m_nvar; k[j++] = 0);
            return k;
        }
        if(i == 0)
            j = deg;
        else {
            for(j = 0; idx_tmp >= m_N[l][j]; ++j);
            k[i-1] -= j;
        }
        k[i] = j;
        idx_tmp -= m_N[l--][j-1];
    }

    return k;
}


template class Chebyshev_Polynomial<double>;
template class Chebyshev_Polynomial<float>;
template class Chebyshev_Polynomial<long double>;
