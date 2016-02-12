/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
---------------- e-mail: annalisa.riccardi@strath.ac.uk ----------------------
------------------------- Author: Annalisa Riccardi --------------------------
*/


#include "Polynomial/base_polynomial.h"

using namespace smart;
using namespace polynomial;

template <class T>
base_polynomial<T>::base_polynomial(const string &name, const int &nvar, const int &order): m_name(""), m_coeffs(0), m_degree(0), m_nvar(0){
    //allocate memory for coefficients vector
    // if(order > MAX_DEGREE){
    //     std::cout<<"Maximum allowed polynomial degree is 20";
    //     exit(EXIT_FAILURE);
    // }

    int n = combination(nvar,order);
    m_coeffs.resize(n);

    //save some info
    m_degree = order;
    m_nvar = nvar;

    m_name = name;

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
base_polynomial<T>::base_polynomial(const string &name, const int &nvar, const int &order, const int &i): m_name(""), m_coeffs(0), m_degree(0), m_nvar(0){

    // if(order > MAX_DEGREE){
    //     std::cout<<"Maximum allowed polynomial degree is 20";
    //     exit(EXIT_FAILURE);
    // }

    //save some info
    m_degree = order;
    m_nvar = nvar;

    m_name = name;

    if (order == 0){
        std::cout<<"cannot assign variables to a polynomial of order 0";
        exit(EXIT_FAILURE);
    }
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
        for(int i=0; i<=nvar; i++){
            m_J[i].resize(order+1);
            m_N[i].resize(order+1);

        }

        initialize_J();
        initialize_N();

    }

}

template <class T>
base_polynomial<T>::base_polynomial(const string &name, const int &nvar, const int &order, const T &value): m_name(""), m_coeffs(0), m_degree(0), m_nvar(0){

    // if(order > MAX_DEGREE){
    //     std::cout<<"Maximum allowed polynomial degree is 20";
    //     exit(EXIT_FAILURE);
    // }

    //save some info
    m_degree = order;
    m_nvar = nvar;

    m_name = name;

    //allocate memory for coefficients vector

    int n = combination(nvar,order);
    m_coeffs.resize(n);
    m_coeffs[0] = value;

    m_J.resize(nvar+1);
    m_N.resize(nvar+1);
    for(int i=0; i<=nvar; i++){
        m_J[i].resize(order+1);
        m_N[i].resize(order+1);

    }

    initialize_J();
    initialize_N();

}

//not part of the algebra, private routines
template <class T>
void base_polynomial<T>::initialize_J()
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
void base_polynomial<T>::initialize_N()
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
int base_polynomial<T>::get_idx(const std::vector<int> &k) const
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
std::vector<int> base_polynomial<T>::get_row(const int &idx, const int &deg) const
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

template class base_polynomial<double>;
template class base_polynomial<float>;
template class base_polynomial<long double>;
