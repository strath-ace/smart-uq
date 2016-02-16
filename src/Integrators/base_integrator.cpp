#include "Integrators/base_integrator.h"

using namespace smart;
using namespace integrator;

template < class T >
base_integrator<T>::base_integrator(const std::string &name, const dynamics::base_dynamics<T> *dyn) : m_name(name), m_dyn(dyn){
}

template < class T >
std::string base_integrator<T>::get_name() const{
    return m_name;
}


template class base_integrator<double>;
template class base_integrator<float>;
template class base_integrator<long double>;
template class base_integrator<polynomial::base_polynomial<double> >;
template class base_integrator<polynomial::base_polynomial<float> >;
template class base_integrator<polynomial::base_polynomial<long double> >;

