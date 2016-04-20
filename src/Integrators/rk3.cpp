/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "../../include/Integrators/rk3.h"


using namespace smartuq;
using namespace smartuq::integrator;

template < class T >
rk3<T>::rk3(const dynamics::base_dynamics<T> *dyn) : base_integrator<T>("Runge Kutta 3 fixed step time", dyn)
{
}

template < class T >
rk3<T>::~rk3(){

}

template < class T >
int rk3<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const{

        // sanity checks
        if(ti<0 || tend<0)
                smart_throw(m_name+": initial time and final time must be greater or equal to 0");
        if(tend<ti)
                smart_throw(m_name+": final time must be greater than initial time");

        xfinal.clear();

        std::vector<T> x(x0), xtemp(x0), k1, k2, k3;

	unsigned int n = x0.size();
	double h = (tend-ti)/nsteps;

    for(int i=0; i<nsteps; i++){
		double t1, t2, t3;
		t1 = ti + i*h;
		t2 = t1 + h/2.0;
		t3 = t1 + h*3.0/4.0;

		//* Evaluate k1 = f(x).
		m_dyn->evaluate(t1, x, k1);

		//* Evaluate k2 = f(x+h/2*k1),
		for(int j=0; j<n; j++)
		    xtemp[j] = x[j]+k1[j]*h/2.0;
		m_dyn->evaluate(t2, xtemp, k2);

		//* Evaluate k3 = f(x+3/4*k2),
		for(int j=0; j<n; j++)
		    xtemp[j] = x[j]+k2[j]*h*3.0/4.0;
		m_dyn->evaluate(t3, xtemp, k3);

		//* Return x(t+h) computed from third-order Runge Kutta.
		for(int j=0; j<n; j++)
		    x[j] += (2.0*k1[j]+3.0*k2[j]+4.0*k3[j])*h/9.0;

	}

	for(int i=0; i<x0.size(); i++)
	    xfinal.push_back(x[i]);

	return 0;
}


template class rk3<double>;
template class rk3<float>;
template class rk3<long double>;
template class rk3<polynomial::chebyshev_polynomial<double> >;
template class rk3<polynomial::chebyshev_polynomial<float> >;
template class rk3<polynomial::chebyshev_polynomial<long double> >;
template class rk3<polynomial::taylor_polynomial<double> >;
template class rk3<polynomial::taylor_polynomial<float> >;
template class rk3<polynomial::taylor_polynomial<long double> >;

