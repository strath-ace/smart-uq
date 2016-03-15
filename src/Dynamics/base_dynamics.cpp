/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "Dynamics/base_dynamics.h"

using namespace smart;
using namespace dynamics;

template < class T >
base_dynamics<T>::base_dynamics(const std::string &name): m_name(name){
}

template < class T >
std::string base_dynamics<T>::get_name() const{
    return m_name;
}

template < class T >
base_dynamics<T>::~base_dynamics(){

}

template class base_dynamics<double>;
template class base_dynamics<float>;
template class base_dynamics<long double>;
template class base_dynamics<polynomial::chebyshev_polynomial<double> >;
template class base_dynamics<polynomial::chebyshev_polynomial<float> >;
template class base_dynamics<polynomial::chebyshev_polynomial<long double> >;
template class base_dynamics<polynomial::taylor_polynomial<double> >;
template class base_dynamics<polynomial::taylor_polynomial<float> >;
template class base_dynamics<polynomial::taylor_polynomial<long double> >;
