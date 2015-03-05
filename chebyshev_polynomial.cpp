/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
---------------- e-mail: annalisa.riccardi@strath.ac.uk ----------------------
------------------------- Author: Annalisa Riccardi --------------------------
*/


#include "chebyshev_polynomial.h"

template <class T>
Chebyshev_Polynomial<T>::Chebyshev_Polynomial(int nvar, int order): m_coeffs(0), m_degree(0), m_nvar(0){
    //allocate memory for coefficients vector
    m_coeffs.resize(factorial(nvar+order)/(factorial(nvar)*factorial(order)));

    //save some info
    m_degree = order;
    m_nvar = nvar;
}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator+(const Chebyshev_Polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(m_degree!=other.get_degree()()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    int n = get_ncoeffs();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    Chebyshev_Polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs + other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Chebyshev_Polynomial<T> Chebyshev_Polynomial<T>::operator-(const Chebyshev_Polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(m_degree!=other.get_degree()()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    int n = get_ncoeffs();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    Chebyshev_Polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs - other_coeffs[i];

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
