/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "../../include/Integrators/base_integrator.h"

using namespace smartuq;
using namespace smartuq::integrator;

template < class T >
base_integrator<T>::base_integrator(const std::string &name, const dynamics::base_dynamics<T> *dyn) : m_name(name), m_dyn(dyn){
}

template < class T >
std::string base_integrator<T>::get_name() const{
    return m_name;
}

template < class T >
base_integrator<T>::~base_integrator(){

}


template class base_integrator<double>;
template class base_integrator<float>;
template class base_integrator<long double>;
template class base_integrator<polynomial::chebyshev_polynomial<double> >;
template class base_integrator<polynomial::chebyshev_polynomial<float> >;
template class base_integrator<polynomial::chebyshev_polynomial<long double> >;
template class base_integrator<polynomial::taylor_polynomial<double> >;
template class base_integrator<polynomial::taylor_polynomial<float> >;
template class base_integrator<polynomial::taylor_polynomial<long double> >;
