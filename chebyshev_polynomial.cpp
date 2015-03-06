/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
---------------- e-mail: annalisa.riccardi@strath.ac.uk ----------------------
------------------------- Author: Annalisa Riccardi --------------------------
*/


#include "chebyshev_polynomial.h"

template <class T>
Chebyshev_Polynomial<T>::Chebyshev_Polynomial(int nvar, int order): m_coeffs(0), m_degree(0), m_nvar(0){
    //allocate memory for coefficients vector

    int n = factorial(nvar+order)/(factorial(nvar)*factorial(order));
    m_coeffs.resize(n);

    //save some info
    m_degree = order;
    m_nvar = nvar;

    m_J.resize(nvar+1);
    m_N.resize(nvar+1);
    for(int i=0; i<=nvar; i++){
        m_J[i].resize(order+1);
        m_N[i].resize(order+1);
    }

    initialize_J();
    initialize_N();

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

}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator/(const Chebyshev_Polynomial<T> &other) const{

}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator+(const T& other) const{

    std::vector<T> coeffs(get_ncoeffs());
    Chebyshev_Polynomial<T> res(m_nvar,m_degree);

    coeffs[0] += other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator-(const T& other) const{

    std::vector<T> coeffs(get_ncoeffs());
    Chebyshev_Polynomial<T> res(m_nvar,m_degree);

    coeffs[0] += other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator*(const T& other) const{

    int n = get_ncoeffs();

    std::vector<T> coeffs(n);
    Chebyshev_Polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] *= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator/(const T& other) const{

    int n = get_ncoeffs();

    std::vector<T> coeffs(n);
    Chebyshev_Polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
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

//template <class T>
//void Chebyshev_Polynomial<T>::initialize_t(){
//    const unsigned int n = std::pow( 2.0, m_nvar);

//}

template <class T>
int Chebyshev_Polynomial<T>::get_idx(const std::vector<int> &k, const int &deg) const
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
