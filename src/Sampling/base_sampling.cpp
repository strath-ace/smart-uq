/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "../../include/Sampling/base_sampling.h"

using namespace smartuq;
using namespace sampling;

template <class T>
base_sampling<T>::base_sampling(const unsigned int &dim, const std::vector<T>& a, const std::vector<T>& b, const string &name): m_name(name), m_dim(dim), m_a(a), m_b(b){
    if(m_a.size() != m_dim)
        smart_throw("Base sampling: hypercube lower bound need to have the same size of the hypercube dimension");
    if(m_b.size() != m_dim)
        smart_throw("Base sampling: hypercube upper bound need to have the same size of the hypercube dimension");
    for(unsigned int i=0; i<m_dim; i++){
        if(m_a[i]>m_b[i])
            smart_throw("Base sampling: hypercube bounds need to be a<b");
    }
}

template <class T>
base_sampling<T>::~base_sampling(){

}

template <class T>
std::vector<T> base_sampling<T>::map(const std::vector<T> &point) const{
    if(point.size()!=m_dim)
        smart_throw("Base sampling: mapping of a point of wrong dimension");

    std::vector<T> res(m_dim);
    for(unsigned int i=0; i<m_dim;i++){
        res[i] = point[i]*(m_b[i]-m_a[i]) + m_a[i];
    }
    return res;
}

template class base_sampling<double>;
template class base_sampling<float>;
template class base_sampling<long double>;
