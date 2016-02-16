#include "Integrators/euler.h"


using namespace smart;
using namespace integrator;

template < class T >
euler<T>::euler(const dynamics::base_dynamics<T> *dyn) : base_integrator<T>("Euler Explicit integration scheme", dyn)
{
}

template < class T >
euler<T>::~euler(){

}

template < class T >
int euler<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> x0, std::vector<T> xfinal) const{

	// sanity checks
    	if(ti<0 || tend<0)
        	smart_exception(m_name+": initial time and final time must be greater or equal to 0");
    	if(tend<ti)
        	smart_exception(m_name+": final time must be greater than initial time");   
	if(x0.size() != xfinal.size())
       		smart_exception(m_name+": initial state must have the same size as final state");   

	std::vector<double> dx(x0.size());
	std::vector<double> x(x0.size());

	double h = (tend-ti)/nsteps;
	x = x0;
	for(int i=0; i<nsteps+1; i++){
		m_dyn->evaluate(ti+i*h, x, dx);
		for(size_t j=0; j<x.size(); j++){
			x[j] += h*dx[j];
		}
	}
	xfinal = x;
	return 0;
}

