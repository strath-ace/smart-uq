/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "Dynamics/vanderpol.h"


using namespace smart;
using namespace dynamics;

template < class T >
vanderpol<T>::vanderpol(const T &mu) : base_dynamics<T>("Van der Pol dynamical system"), m_mu(mu)
{

}

template < class T >
vanderpol<T>::~vanderpol()
{

}

template < class T >
int vanderpol<T>::evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const
{
    //sanity checks
    if(t<0)
        smart_throw(m_name+": negative time supplied in evaluation of the dynamical system");
    if(state.size()!=2)
        smart_throw(m_name+": the state dimension needs to be 2");

    dstate.clear();

    dstate.push_back(state[1]);
    dstate.push_back(m_mu*(1.0-state[0]*state[0])*state[1]-state[0]);

    return 0;
}

template class vanderpol<double>;
template class vanderpol<float>;
template class vanderpol<long double>;
template class vanderpol<polynomial::chebyshev_polynomial<double> >;
template class vanderpol<polynomial::chebyshev_polynomial<float> >;
template class vanderpol<polynomial::chebyshev_polynomial<long double> >;
template class vanderpol<polynomial::taylor_polynomial<double> >;
template class vanderpol<polynomial::taylor_polynomial<float> >;
template class vanderpol<polynomial::taylor_polynomial<long double> >;


