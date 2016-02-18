#include "Integrators/rk4.h"


using namespace smart;
using namespace integrator;

template < class T >
rk4<T>::rk4(const dynamics::base_dynamics<T> *dyn) : base_integrator<T>("Runge Kutta 4 fixed step time", dyn)
{
}

template < class T >
rk4<T>::~rk4(){

}

template < class T >
int rk4<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> x0, std::vector<T> xfinal) const{

        // sanity checks
        if(ti<0 || tend<0)
                smart_exception(m_name+": initial time and final time must be greater or equal to 0");
        if(tend<ti)
                smart_exception(m_name+": final time must be greater than initial time");
        if(x0.size() != xfinal.size())
                smart_exception(m_name+": initial state must have the same size as final state");

        std::vector<double> x(x0.size()), xtemp(x0.size()), k1(x0.size()), k2(x0.size()), k3(x0.size()), k4(x0.size());

	unsigned int n = x0.size();
	double h = (tend-ti)/nsteps;
	x = x0;

	for(int i=0; i<nsteps+1; i++){
		//* Evaluate k1 = f(x).
		m_dyn->evaluate(ti+i*h, x, k1);

		//* Evaluate k2 = f(x+h/2*k1),
		for(int j=0; j<n; j++)
		    xtemp[j] = x[j]+k1[j]*h/2.0;
		m_dyn->evaluate(ti+i*h, xtemp, k2);

		//* Evaluate k3 = f(x+h/2*k2),
		for(int j=0; j<n; j++)
		    xtemp[j] = x[j]+k2[j]*h/2.0;
		m_dyn->evaluate(ti+i*h, xtemp, k3);

		//* Evaluate k4 = f(x+h*k3),
		for(int j=0; j<n; j++)
		    xtemp[j] = x[j]+k3[j]*h;
		m_dyn->evaluate(ti+i*h, xtemp, k4);

		//* Return x(t+h) computed from second-order Runge Kutta.
		for(int j=0; j<n; j++)
		    x[j] += (k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j])*h/6.0;

	}

	xfinal = x;

	return 0;
}


