/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
---------------- e-mail: annalisa.riccardi@strath.ac.uk ----------------------
------------------------- Author: Annalisa Riccardi --------------------------
*/


#include "chebyshev_polynomial.h"

template <class T>
Chebyshev_Polynomial<T>::Chebyshev_Polynomial(const int &nvar, const int &order): m_coeffs(0), m_degree(0), m_nvar(0){
    //allocate memory for coefficients vector

    int n = factorial(nvar+order)/(factorial(nvar)*factorial(order));
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

    //save some info
    m_degree = order;
    m_nvar = nvar;

    if(i>=m_nvar){
        std::cout<<"base elements index are from [0,nvar-1]";
        exit(EXIT_FAILURE);
    }
    else{
        //allocate memory for coefficients vector

        int n = factorial(nvar+order)/(factorial(nvar)*factorial(order));
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

    //save some info
    m_degree = order;
    m_nvar = nvar;

        //allocate memory for coefficients vector

        int n = factorial(nvar+order)/(factorial(nvar)*factorial(order));
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
    std::vector<T> res_coeffs(m_coeffs.size());
    double nvariations = pow(2,m_nvar);

    std::vector<T> other_coeffs = other.get_coeffs();
    for(int i=0; i<=m_degree; i++){//loop over subset degree i of poly1
        for(int j=0; j<=other.get_degree(); j++){//loop over subset degree j of poly2
              if((i+j)<=m_degree){
                  for(int idx1=0; idx1<m_J[m_nvar][i]; idx1++){//index over elements with degree i in poly1
                      for(int idx2=0; idx2<m_J[m_nvar][j]; idx2++){//index over elements with degree j in poly2
                          std::vector<int> v1 = get_row(idx1,i);
                          std::vector<int> v2 = get_row(idx2,j);
                          std::vector<int> v3(m_nvar);
                          for(int iter=0; iter<nvariations; iter++){
                              for(int k=0; k<m_nvar; k++){
                                  v3[k] = std::fabs(v1[k]+m_t[iter][k]*v2[k]);
                              }
                              int pos = get_idx(v3);
                              int deg3 = std::accumulate(v3.begin(),v3.end(),0);
                              int sub_idx1=0, sub_idx2=0, sub_idx3=0;
                              if(deg3>0) sub_idx1=m_N[m_nvar][deg3-1];
                              if(i>0) sub_idx2=m_N[m_nvar][i-1];
                              if(j>0) sub_idx3=m_N[m_nvar][j-1];

                              res_coeffs[sub_idx1 + pos] +=
                                      (1.0/nvariations)*(m_coeffs[sub_idx2+idx1]*other_coeffs[sub_idx3+idx2]);
                          }
                      }
                  }
              }
        }
    }

    res.set_coeffs(res_coeffs);
    return res;
}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator/(const Chebyshev_Polynomial<T> &other) const{

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
