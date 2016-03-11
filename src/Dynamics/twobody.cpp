/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "Dynamics/twobody.h"


using namespace smart;
using namespace dynamics;

template < class T >
twobody<T>::twobody(const std::vector<T> &param, const double &t_scale, const double &r_scale) : base_dynamics<T>("Two Body Problem"),
    m_param(param), m_t_scale(t_scale), m_r_scale(r_scale)
{
    if(m_param.size()!=10)
        smart_throw(m_name+": the parameters list need to be of size 10");

}

template < class T >
twobody<T>::~twobody()
{

}

template < class T >
int twobody<T>::evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const
{
    //sanity checks
    if(t<0)
        smart_throw(m_name+": negative time supplied in evaluation of the dynamical system");
    if(state.size()!=7)
        smart_throw(m_name+": the state dimension needs to be 7");

    dstate.clear();

    //constant parameters
    double radius_earth = 6378*pow(10,3) / m_r_scale;
    double mu_earth = 398600.4415*pow(10,9) / (pow(m_r_scale,3)/pow(m_t_scale,2));
    double omega_earth = 7.2921150*pow(10,-5) / (1.0/m_t_scale);
    double H0_atmosphere = 900000 / m_r_scale;

    //precomputations
    T r = sqrt(state[0]*state[0]+state[1]*state[1]+state[2]*state[2]);
    T tmp_3D =  mu_earth/pow(r, 3);

    //atmospheric model
    T rho = m_param[4]*exp(-(r-radius_earth-H0_atmosphere)/m_param[5]);

    //relative velocity
    T rel_v_x = state[3]-omega_earth*state[1];
    T rel_v_y = state[4]+omega_earth*state[0];
    T mod_rel_v = sqrt(rel_v_x*rel_v_x+rel_v_y*rel_v_y+state[5]*state[5]);

    //drag computation
    T tmp_drag = 0.5*rho*m_param[6]*mod_rel_v/state[6];

    dstate.push_back(state[3]); //dx/dt
    dstate.push_back(state[4]); //dy/dt
    dstate.push_back(state[5]); //dz/dt
    dstate.push_back(-tmp_3D*state[0]+m_param[0]/state[6]+m_param[7]-tmp_drag*rel_v_x); //dv_x/dt
    dstate.push_back(-tmp_3D*state[1]+m_param[1]/state[6]+m_param[8]-tmp_drag*rel_v_y); //dv_y/dt
    dstate.push_back(-tmp_3D*state[2]+m_param[2]/state[6]+m_param[9]-tmp_drag*state[5]); //dv_z/dt
    dstate.push_back(-m_param[3]*sqrt(m_param[0]*m_param[0]+m_param[1]*m_param[1]+m_param[2]*m_param[2])); //dm/dt

    return 0;

}


template class twobody<double>;
template class twobody<float>;
template class twobody<long double>;
template class twobody<polynomial::chebyshev_polynomial<double> >;
template class twobody<polynomial::chebyshev_polynomial<float> >;
template class twobody<polynomial::chebyshev_polynomial<long double> >;
template class twobody<polynomial::taylor_polynomial<double> >;
template class twobody<polynomial::taylor_polynomial<float> >;
template class twobody<polynomial::taylor_polynomial<long double> >;

